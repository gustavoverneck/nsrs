// Utilitários de linha de comando compartilhados pelos subcomandos.

use std::collections::HashMap;

use nsrs::core::model::ModelParams;
use nsrs::{FSU2, GM1, GM3};

/// Argumentos de um subcomando: posicionais, opções `--nome valor` e
/// chaves booleanas `--nome` (declaradas em `switches`).
pub struct Args {
    pub positional: Vec<String>,
    options: HashMap<String, String>,
    switches: Vec<String>,
}

impl Args {
    pub fn parse(raw: &[String], switches: &[&str]) -> Result<Self, String> {
        let mut args = Args { positional: Vec::new(), options: HashMap::new(), switches: Vec::new() };
        let mut iter = raw.iter();
        while let Some(arg) = iter.next() {
            if let Some(name) = arg.strip_prefix("--") {
                if switches.contains(&name) {
                    args.switches.push(name.to_string());
                } else {
                    let value = iter.next().ok_or_else(|| format!("falta o valor de --{name}"))?;
                    args.options.insert(name.to_string(), value.clone());
                }
            } else {
                args.positional.push(arg.clone());
            }
        }
        Ok(args)
    }

    pub fn switch(&self, name: &str) -> bool {
        self.switches.iter().any(|s| s == name)
    }

    pub fn value(&self, name: &str) -> Option<&str> {
        self.options.get(name).map(String::as_str)
    }

    pub fn f64_or(&self, name: &str, default: f64) -> Result<f64, String> {
        match self.value(name) {
            Some(v) => v.parse().map_err(|_| format!("--{name}: número inválido '{v}'")),
            None => Ok(default),
        }
    }

    pub fn usize_or(&self, name: &str, default: usize) -> Result<usize, String> {
        match self.value(name) {
            Some(v) => v.parse().map_err(|_| format!("--{name}: inteiro inválido '{v}'")),
            None => Ok(default),
        }
    }

    /// Número de threads: `--threads N` ou todos os núcleos disponíveis.
    pub fn threads(&self) -> Result<usize, String> {
        let default = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
        self.usize_or("threads", default)
    }

    /// Modelos de `--models GM1,GM3,...` (padrão: `default`).
    pub fn models(&self, default: &[&str]) -> Result<Vec<(String, ModelParams)>, String> {
        let names: Vec<String> = match self.value("models") {
            Some(list) => list.split(',').map(|s| s.trim().to_string()).collect(),
            None => default.iter().map(|s| s.to_string()).collect(),
        };
        names.into_iter().map(|n| model(&n).map(|m| (n, m))).collect()
    }
}

pub fn model(name: &str) -> Result<ModelParams, String> {
    match name {
        "GM1" => Ok(GM1),
        "GM3" => Ok(GM3),
        "FSU2" => Ok(FSU2),
        _ => Err(format!("modelo desconhecido '{name}' (use GM1, GM3 ou FSU2)")),
    }
}

pub fn f64_list(items: &[String], what: &str) -> Result<Vec<f64>, String> {
    items
        .iter()
        .map(|s| s.parse().map_err(|_| format!("{what}: número inválido '{s}'")))
        .collect()
}

/// Notação científica com duas casas, sem '+': 1.00e17, 1.00e-3.
pub fn format_sci(value: f64) -> String {
    format!("{:.2e}", value).replace('+', "")
}

pub fn create_dir(path: &str) -> Result<(), String> {
    std::fs::create_dir_all(path).map_err(|e| format!("não foi possível criar '{path}': {e}"))
}
