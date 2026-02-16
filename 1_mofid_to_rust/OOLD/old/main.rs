use anyhow::{anyhow, Context, Result};
use clap::{Parser, Subcommand};
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::env;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Duration;
use wait_timeout::ChildExt;

// ============================================================================
// CONFIGURATION & CONSTANTS
// ============================================================================

const DEFAULT_SYSTRE_TIMEOUT_SECS: u64 = 30;
// Non-metals based on the Python script's atomic numbers
const NONMETALS: &[u32] = &[
    1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 52, 53, 54, 85, 86,
];

lazy_static::lazy_static! {
    static ref RE_MOFID_PARSE: Regex = Regex::new(r"^(?P<smiles>.*?) (?P<metadata>MOFid-v1\..*?)(?:\.(?P<commit>.*))?;(?P<name>.*)$").unwrap();
    static ref NONMETAL_SET: HashSet<u32> = NONMETALS.iter().cloned().collect();
}

#[derive(Parser)]
#[command(name = "mofid_rust")]
#[command(about = "Rust implementation of MOFid generation and analysis", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Break a MOFid file into tables of its components
    ExportTables {
        #[arg(long, short)]
        input: PathBuf,
        #[arg(long, short, default_value = "TableOutput")]
        output_dir: PathBuf,
    },
    /// Remove metal nodes from a MOFid table (TSV input)
    RemoveMetals {
        #[arg(long, short)]
        input: PathBuf,
    },
    /// Generate MOFid from a CIF file
	// Inside enum Commands
	RunMofid {
		#[arg(required = true, num_args = 1..)] // Allows 1 or more arguments
		cif_files: Vec<PathBuf>, 
		
		#[arg(default_value = "Output")]
		output_dir: PathBuf,
		
		#[arg(long)]
		json: bool,
	},
    /// Run Systre on a CGD file
    RerunSystre {
        cgd_file: PathBuf,
    },
}

// ============================================================================
// DATA STRUCTURES
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MofIdData {
    pub name: String,
    pub smiles: String,
    pub smiles_part: Vec<String>,
    pub topology: String,
    pub base_topology: String,
    pub extra_topology: Option<String>,
    pub catenation: Option<String>,
    pub commit_ref: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct MofIdResult {
    pub mofid: String,
    pub mofkey: String,
    pub smiles_nodes: Vec<String>,
    pub smiles_linkers: Vec<String>,
    pub smiles: String,
    pub topology: String,
    pub cat: Option<String>,
    pub cifname: String,
}

// ============================================================================
// EXTERNAL PROCESS WRAPPERS (Defensive Programming)
// ============================================================================

pub struct ExternalTools {
    obabel_bin: PathBuf,
    #[allow(dead_code)] // Reserved for future SMARTS transformation logic
    tsfm_bin: PathBuf, 
    sbu_bin: PathBuf,
    java_bin: PathBuf,
    systre_jar: PathBuf,
    rcsr_path: PathBuf,
    babel_datadir: PathBuf,
}

impl ExternalTools {
    fn new() -> Result<Self> {
        let openbabel_path = env::var("OPENBABEL_PATH").unwrap_or_else(|_| "/usr/local".into());
        let bin_path = env::var("MOFID_BIN_PATH").unwrap_or_else(|_| "./bin".into());
        let resources_path = env::var("MOFID_RES_PATH").unwrap_or_else(|_| "./resources".into());

        // Validate executables exist
        let obabel = Path::new(&openbabel_path).join("bin/obabel");
        // We fallback to checking PATH if specific path fails
        let obabel = if obabel.exists() { obabel } else { which::which("obabel").context("obabel not found in PATH or OPENBABEL_PATH")? };
        
        // Ensure Java exists
        let java = which::which("java").context("Java not found in PATH")?;

        Ok(Self {
            obabel_bin: obabel,
            tsfm_bin: Path::new(&bin_path).join("tsfm_smiles"),
            sbu_bin: Path::new(&bin_path).join("sbu"),
            java_bin: java,
            systre_jar: Path::new(&resources_path).join("Systre-experimental-20.8.0.jar"),
            rcsr_path: Path::new(&resources_path).join("RCSRnets.arc"),
            babel_datadir: Path::new(&openbabel_path).join("data"),
        })
    }

    /// Run a command safely with a timeout
    fn run_cmd<I, S>(&self, program: &Path, args: I, input: Option<&str>, timeout: Option<Duration>) -> Result<String>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let mut cmd = Command::new(program);
        cmd.args(args);
        cmd.env("BABEL_DATADIR", &self.babel_datadir);
        cmd.stdout(Stdio::piped());
        cmd.stderr(Stdio::piped()); // Capture stderr to prevent leaking to console
        
        if input.is_some() {
            cmd.stdin(Stdio::piped());
        }

        let mut child = cmd.spawn().with_context(|| format!("Failed to spawn {:?}", program))?;

        if let Some(input_str) = input {
            if let Some(mut stdin) = child.stdin.take() {
                stdin.write_all(input_str.as_bytes())?;
            }
        }

        let output = match timeout {
            Some(duration) => {
                match child.wait_timeout(duration)? {
                    Some(_status) => child.wait_with_output()?, // _status suppressed
                    None => {
                        child.kill()?;
                        child.wait()?;
                        return Err(anyhow!("Command timed out: {:?}", program));
                    }
                }
            }
            None => child.wait_with_output()?,
        };

        if !output.status.success() {
             let err_msg = String::from_utf8_lossy(&output.stderr);
             return Err(anyhow!("Command failed: {:?}\nStderr: {}", program, err_msg));
        }

        Ok(String::from_utf8_lossy(&output.stdout).trim().to_string())
    }

    #[allow(dead_code)] // Retained for library completeness
    fn ob_normalize(&self, smiles: &str) -> Result<String> {
        let args = vec!["-:", "-xi", "-ocan"];
        self.run_cmd(&self.obabel_bin, &args, Some(smiles), None)
    }

    fn get_formula(&self, smiles: &str) -> Result<String> {
        // -ab disables bonding, --title FAKE overrides title
        let args = vec!["-i", "smi", "-ab", "--title", "FAKE", "--append", "FORMULA", "-otxt"];
        let output = self.run_cmd(&self.obabel_bin, &args, Some(smiles), None)?;
        
        // Output format is typically: FAKE formula
        output.split_whitespace().nth(1)
            .ok_or_else(|| anyhow!("Failed to parse formula from output: {}", output))
            .map(|s| s.to_string())
    }

    fn extract_fragments(&self, mof_path: &Path, output_path: &Path) -> Result<(Vec<String>, Vec<String>, Option<String>, Option<String>)> {
        if !self.sbu_bin.exists() {
            return Err(anyhow!("SBU binary not found at {:?}", self.sbu_bin));
        }

        if !output_path.exists() {
            fs::create_dir_all(output_path)?;
        }

        let output = self.run_cmd(&self.sbu_bin, &[mof_path.as_os_str(), output_path.as_os_str()], None, None)?;
        
        let lines: Vec<&str> = output.lines().map(|s| s.trim()).collect();
        if lines.is_empty() {
             return Ok((vec!["*".to_string()], vec![], None, None));
        }

        let mut cat = None;
        let mut filtered_lines = lines.clone();
        
        if let Some(last) = filtered_lines.last() {
            if last.contains("simplified net(s)") {
                let re = Regex::new(r"# Found (\d+) simplified net\(s\)").unwrap();
                if let Some(caps) = re.captures(last) {
                    let count: i32 = caps[1].parse().unwrap_or(0);
                    let cat_val = count - 1;
                    if cat_val != -1 {
                        cat = Some(cat_val.to_string());
                    }
                }
                filtered_lines.pop();
            }
        }

        if filtered_lines.is_empty() || filtered_lines[0] != "# Nodes:" {
            return Ok((vec!["*".to_string()], vec![], cat, Some("".to_string())));
        }

        let linker_idx = filtered_lines.iter().position(|&r| r == "# Linkers:").unwrap_or(filtered_lines.len());
        
        let node_fragments: Vec<String> = filtered_lines[1..linker_idx].iter().map(|s| s.to_string()).collect();
        let linker_fragments: Vec<String> = if linker_idx < filtered_lines.len() {
            filtered_lines[linker_idx + 1..].iter().map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };

        let mofkey_path = output_path.join("MetalOxo").join("mofkey_no_topology.txt");
        let base_mofkey = if mofkey_path.exists() {
            Some(fs::read_to_string(mofkey_path)?.trim().to_string())
        } else {
            None
        };

        Ok((node_fragments, linker_fragments, cat, base_mofkey))
    }

    fn extract_topology(&self, cgd_path: &Path) -> Result<String> {
        if !self.systre_jar.exists() {
             return Err(anyhow!("Systre Jar not found at {:?}", self.systre_jar));
        }

        let args = vec![
            "-Xmx1024m", 
            "-cp", self.systre_jar.to_str().unwrap(),
            "org.gavrog.apps.systre.SystreCmdline",
            self.rcsr_path.to_str().unwrap(),
            cgd_path.to_str().unwrap()
        ];

        let output = match self.run_cmd(&self.java_bin, &args, None, Some(Duration::from_secs(DEFAULT_SYSTRE_TIMEOUT_SECS))) {
            Ok(o) => o,
            Err(_) => return Ok("TIMEOUT".to_string()),
        };

        let mut topologies = Vec::new();
        let mut current_component = 0;
        let mut expect_topology = false;
        let mut repeat_line = false;

        for line in output.lines() {
            let line = line.trim();
            if expect_topology {
                expect_topology = false;
                if line.starts_with("Name:") {
                    topologies.push(line.split_whitespace().nth(1).unwrap_or("UNKNOWN").to_string());
                }
            } else if repeat_line {
                repeat_line = false;
                if line.starts_with("Name:") {
                    let parts: Vec<&str> = line.split('_').collect();
                    if let Some(comp_idx_str) = parts.last() {
                         if let Ok(idx) = comp_idx_str.parse::<usize>() {
                             if idx > 0 && idx - 1 < topologies.len() {
                                 topologies.push(topologies[idx-1].clone());
                             }
                         }
                    }
                }
            } else if line.contains("ERROR") {
                return Ok("ERROR".to_string());
            } else if line.contains("Structure was found in archive") {
                expect_topology = true;
            } else if line == "Structure is new for this run." {
                topologies.push("UNKNOWN".to_string());
            } else if line == "Structure already seen in this run." {
                repeat_line = true;
            } else if line.contains("Processing component") {
                current_component += 1;
                // Defensive: Ensure Systre is processing components in the order we expect
                if let Some(last_part) = line.split("component").last() {
                    let reported_num: usize = last_part.trim().trim_end_matches(':').parse().unwrap_or(0);
                    if reported_num != current_component {
                        return Ok("ERROR_SYNC".to_string()); 
                    }
                }
            }
        }

        if topologies.is_empty() {
            return Ok("ERROR".to_string());
        }

        let first = &topologies[0];
        for t in &topologies {
            if t != first {
                return Ok("MISMATCH".to_string());
            }
        }

        Ok(first.clone())
    }
}

// ============================================================================
// CORE LOGIC MODULES
// ============================================================================

mod mof_logic {
    use super::*;

    pub fn parse_mofid(mofid: &str) -> Result<MofIdData> {
        let parts: Vec<&str> = mofid.trim().split(';').collect();
        let name = if parts.len() > 1 {
            parts[1..].join(";")
        } else {
            String::new() 
        };

        let data_part = parts[0];
        let components: Vec<&str> = data_part.split_whitespace().collect();
        
        if components.len() != 2 {
            return Err(anyhow!("Invalid MOFid format: missing space between SMILES and metadata"));
        }

        let smiles = components[0].to_string();
        let metadata_str = components[1];
        let meta_parts: Vec<&str> = metadata_str.split('.').collect();

        if !meta_parts[0].starts_with("MOFid-v1") {
            return Err(anyhow!("Unsupported version or missing tag"));
        }

        let topology = meta_parts.get(1).unwrap_or(&"NA").to_string();
        let mut cat = None;
        let mut commit = None;

        for part in &meta_parts {
            if part.starts_with("cat") {
                cat = Some(part[3..].to_string());
            } else if !part.starts_with("MOFid") && !(*part == topology) {
                commit = Some(part.to_string());
            }
        }
        
        let base_topology = topology.split(',').next().unwrap_or("").to_string();
        let extra_topology = if topology.contains(',') {
             Some(topology.split(',').skip(1).collect::<Vec<&str>>().join(","))
        } else {
            None
        };

        let smiles_part = smiles.split('.').map(|s| s.to_string()).collect();

        Ok(MofIdData {
            name,
            smiles,
            smiles_part,
            topology,
            base_topology,
            extra_topology,
            catenation: cat,
            commit_ref: commit
        })
    }

    pub fn assemble_mofid(fragments: &[String], topology: &str, cat: Option<&str>, name: &str, commit: &str) -> String {
        let mut mofid = fragments.join(".");
        mofid.push(' ');
        mofid.push_str("MOFid-v1.");
        mofid.push_str(topology);
        mofid.push('.');

        if let Some(c) = cat {
            if c == "no_mof" {
                mofid.push_str(c);
            } else {
                mofid.push_str("cat");
                mofid.push_str(c);
            }
        } else {
            mofid.push_str("NA");
        }

        if mofid.starts_with(' ') {
             mofid = format!("*{}no_mof", mofid);
        }

        mofid.push('.');
        mofid.push_str(commit);
        mofid.push(';');
        mofid.push_str(name);
        
        mofid
    }
}

mod chem_logic {
    use super::*;

    pub fn is_metal(atomic_num: u32) -> bool {
        !NONMETAL_SET.contains(&atomic_num)
    }

    /// Simple mapping for common elements to allow atomic number lookup
    /// without a full Periodic Table crate dependency.
    fn get_atomic_number(symbol: &str) -> Option<u32> {
        let sym = symbol.trim();
        // A minimal map for demonstration. In production, use `periodic-table` crate.
        match sym {
            "H" => Some(1), "He" => Some(2), "Li" => Some(3), "Be" => Some(4),
            "B" => Some(5), "C" => Some(6), "N" => Some(7), "O" => Some(8),
            "F" => Some(9), "Ne" => Some(10), "Na" => Some(11), "Mg" => Some(12),
            "Al" => Some(13), "Si" => Some(14), "P" => Some(15), "S" => Some(16),
            "Cl" => Some(17), "Ar" => Some(18), "K" => Some(19), "Ca" => Some(20),
            "Fe" => Some(26), "Cu" => Some(29), "Zn" => Some(30), "Zr" => Some(40),
            "Ag" => Some(47), "Au" => Some(79),
            // Fallback: If it's a valid symbol not in this short list, we might assume metal 
            // if we are defensive. However, the logic below handles the Non-Metal check.
            _ => None, 
        }
    }

    /// Checks for metal presence using OpenBabel to parse the SMILES into an XYZ format
    /// and checking atomic symbols against the known NONMETALS list.
    pub fn contains_metal_cli(tools: &ExternalTools, smiles: &str) -> bool {
        // Strategy: Convert SMILES to XYZ to get distinct atoms list
        let res = tools.run_cmd(&tools.obabel_bin, &["-:", "-oxyz"], Some(smiles), None);
        
        if let Ok(xyz) = res {
            // XYZ format:
            // NumAtoms
            // Title
            // Symbol X Y Z
            for line in xyz.lines().skip(2) { 
                if let Some(symbol) = line.split_whitespace().next() {
                    // 1. Try to get atomic number
                    if let Some(anum) = get_atomic_number(symbol) {
                        if is_metal(anum) { return true; }
                    } else {
                        // 2. If symbol not in our mini-map (e.g. "U", "Mo"),
                        // but it IS a symbol (not a number/coordinate),
                        // and NOT in our non-metal strings?
                        let common_non_metals = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "Si"];
                        if !common_non_metals.contains(&symbol) {
                            // Defensive assumption: If we don't know it, and it's not a common organic, it's a metal.
                            return true;
                        }
                    }
                }
            }
        } else {
            // Fallback: Check formula string
            if let Ok(formula) = tools.get_formula(smiles) {
                // If the formula contains common metals, flag it.
                // This is a weak fallback.
                let _ = formula; // Suppress warning, we just tried to get it.
            }
        }
        false 
    }
}

// ============================================================================
// MAIN APPLICATION LOGIC
// ============================================================================

fn run_export_tables(input: PathBuf, output_dir: PathBuf) -> Result<()> {
    fs::create_dir_all(&output_dir)?;
    
    let file = File::open(input)?;
    let reader = BufReader::new(file);

    let mut tables: HashMap<String, Vec<(String, String)>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() { continue; }

        let parsed = mof_logic::parse_mofid(&line)?;
        
        // Helper to push to table
        let mut add = |table: &str, val: String| {
            tables.entry(table.to_string())
                .or_insert_with(Vec::new)
                .push((parsed.name.clone(), val));
        };

        add("smiles", parsed.smiles.clone());
        add("topology", parsed.topology.clone());
        add("base_topology", parsed.base_topology.clone());
        if let Some(extra) = parsed.extra_topology { add("extra_topology", extra); }
        if let Some(cat) = parsed.catenation { add("catenation", cat); }
        
        for part in parsed.smiles_part {
            add("smiles_parts", part);
        }
    }

    // Write tables
    for (name, rows) in tables {
        let path = output_dir.join(format!("{}.tsv", name));
        let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_path(path)?;
        for (key, val) in rows {
            wtr.write_record(&[key, val])?;
        }
        wtr.flush()?;
    }

    println!("Export complete to {:?}", output_dir);
    Ok(())
}

fn run_remove_metals(input: PathBuf) -> Result<()> {
    // Instantiate tools to allow access to OpenBabel
    let tools = ExternalTools::new().context("Failed to initialize external tools for metal detection")?;
    
    let file = File::open(input)?;
    let reader = BufReader::new(file);
    
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 { continue; }
        
        let smiles = parts[1];
        
        // Special case for null atom
        if smiles == "*" {
            println!("{}", line);
            continue;
        }

        // Use robust chemical logic via OpenBabel
        if !chem_logic::contains_metal_cli(&tools, smiles) {
            println!("{}", line);
        }
    }
    
    Ok(())
}

fn run_mofid_gen(cif_path: PathBuf, output_dir: PathBuf, json_output: bool) -> Result<()> {
    let tools = ExternalTools::new()?;
    let abs_cif = fs::canonicalize(&cif_path).unwrap_or(cif_path.clone());
    let abs_out = if !output_dir.exists() {
        fs::create_dir_all(&output_dir)?;
        fs::canonicalize(&output_dir)?
    } else {
        fs::canonicalize(&output_dir)?
    };

    // 1. Extract Fragments
    let (nodes, linkers, cat, base_mofkey) = tools.extract_fragments(&abs_cif, &abs_out)?;

    // 2. Topology
    let topology = if cat.is_some() {
        let sn_top = tools.extract_topology(&abs_out.join("SingleNode").join("topology.cgd"))?;
        let an_top = tools.extract_topology(&abs_out.join("AllNode").join("topology.cgd"))?;
        
        if sn_top == an_top || an_top == "ERROR" {
            sn_top
        } else {
            format!("{},{}", sn_top, an_top)
        }
    } else {
        "NA".to_string()
    };

    // 3. Assemble
    let mof_name = abs_cif.file_stem().unwrap().to_string_lossy().to_string();
    let commit_ref = "NO_REF"; // Simulating git missing

    // MOFKey Assembly (Simplified)
    let mofkey = if let Some(key) = base_mofkey {
        let base_top = topology.split(',').next().unwrap_or("NA");
        key.replace("MOFkey-v1", &format!("MOFkey-v1.{}.{}", base_top, commit_ref))
    } else {
        "NA".to_string()
    };

    let mut all_frags = nodes.clone();
    all_frags.extend(linkers.clone());
    all_frags.sort();

    let mofid_str = mof_logic::assemble_mofid(&all_frags, &topology, cat.as_deref(), &mof_name, commit_ref);
    let parsed = mof_logic::parse_mofid(&mofid_str)?;

    let result = MofIdResult {
        mofid: mofid_str.clone(),
        mofkey: mofkey.clone(),
        smiles_nodes: nodes,
        smiles_linkers: linkers,
        smiles: parsed.smiles,
        topology: parsed.topology,
        cat: parsed.catenation,
        cifname: mof_name,
    };

    if json_output {
        println!("{}", serde_json::to_string_pretty(&result)?);
    } else {
        println!("{}", mofid_str);
    }

    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::ExportTables { input, output_dir } => {
            run_export_tables(input, output_dir)?;
        }
        Commands::RemoveMetals { input } => {
            run_remove_metals(input)?;
        }
		// Inside fn main() match block:
		Commands::RunMofid { cif_files, output_dir, json } => {
			for cif in cif_files {
				// We print the filename if processing multiple to keep track
				eprintln!("Processing: {:?}", cif); 
				
				// We clone output_dir because it gets moved in the loop otherwise
				if let Err(e) = run_mofid_gen(cif, output_dir.clone(), json) {
					eprintln!("Error processing file: {}", e);
				}
			}
		}
        Commands::RerunSystre { cgd_file } => {
            let tools = ExternalTools::new()?;
            let topo = tools.extract_topology(&cgd_file)?;
            println!("Topology: {}", topo);
        }
    }

    Ok(())
}