use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::Duration;
use std::env;
use std::ffi::OsStr;
use anyhow::{Result, anyhow, Context};
use wait_timeout::ChildExt;

/// Manages external binaries required for MOF decomposition.
pub struct ExternalTools {
    pub obabel_bin: PathBuf,
    pub sbu_bin: PathBuf,
    pub java_bin: PathBuf,
    pub systre_jar: PathBuf,
    pub babel_datadir: PathBuf,
}

impl ExternalTools {
    /// Locates all necessary binaries in the system path or environment variables.
    pub fn new() -> Result<Self> {
        let current_dir = env::current_dir().context("Cannot determine current directory")?;
        
        // Environment overrides -> Local defaults -> System PATH
        let bin_path = env::var("MOFID_BIN_PATH").map(PathBuf::from).unwrap_or_else(|_| current_dir.join("bin"));
        let res_path = env::var("MOFID_RES_PATH").map(PathBuf::from).unwrap_or_else(|_| current_dir.join("resources"));
        let obabel_path = env::var("OPENBABEL_PATH").map(PathBuf::from).unwrap_or_else(|_| current_dir.clone());

        // Locate OpenBabel
        let obabel = which::which("obabel").or_else(|_| {
            // Fallback for Docker/Linux standard paths
            if Path::new("/usr/bin/obabel").exists() { Ok(PathBuf::from("/usr/bin/obabel")) }
            else { Err(anyhow!("'obabel' not found in PATH.")) }
        })?;
        
        // Locate Java
        let java = which::which("java").context("Java not found in PATH")?;

        let tools = Self {
            obabel_bin: obabel,
            sbu_bin: bin_path.join("sbu"), // The python/bash script that splits the MOF
            java_bin: java,
            systre_jar: res_path.join("Systre-experimental-20.8.0.jar"),
            babel_datadir: obabel_path.join("data"),
        };

        // Validation
        if !tools.sbu_bin.exists() { 
            return Err(anyhow!("SBU binary missing at {:?}. Please set MOFID_BIN_PATH.", tools.sbu_bin)); 
        }
        if !tools.systre_jar.exists() { 
            // Warning only? No, strictly required for topology.
            // But for simple splitting, maybe not? Let's verify strictness.
            // Return error for now to be safe.
             return Err(anyhow!("Systre JAR missing at {:?}. Please set MOFID_RES_PATH.", tools.systre_jar)); 
        }

        Ok(tools)
    }

    /// Executes a command with timeout and environment setup.
    pub fn run_cmd<I, S>(&self, program: &Path, args: I, input: Option<&str>, timeout: Option<Duration>) -> Result<String>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let mut cmd = Command::new(program);
        cmd.args(args);
        
        // Essential OpenBabel env vars
        cmd.env("BABEL_DATADIR", &self.babel_datadir);
        if let Ok(libdir) = env::var("BABEL_LIBDIR") {
            cmd.env("BABEL_LIBDIR", libdir);
        }

        cmd.stdout(Stdio::piped());
        cmd.stderr(Stdio::piped()); 
        
        if input.is_some() {
            cmd.stdin(Stdio::piped());
        }

        let mut child = cmd.spawn().with_context(|| format!("Failed to spawn {:?}", program))?;

        if let Some(input_str) = input {
            if let Some(mut stdin) = child.stdin.take() {
                // Ignore write errors (broken pipe if process exits early)
                let _ = std::io::Write::write_all(&mut stdin, input_str.as_bytes());
            }
        }

        let output = match timeout {
            Some(duration) => {
                match child.wait_timeout(duration)? {
                    Some(_status) => child.wait_with_output()?,
                    None => {
                        let _ = child.kill();
                        let _ = child.wait();
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
}