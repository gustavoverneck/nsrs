// Utilitários de linha de comando compartilhados pelos subcomandos.

use std::collections::HashMap;
use std::time::Duration;

use indicatif::{ProgressBar, ProgressStyle};

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

/// Frases da barra de progresso dos estudos; trocam a cada poucos segundos.
const PROGRESS_PHRASES: [&str; 12] = [
    "colapsando estrelas",
    "empilhando níveis de Landau",
    "convencendo nêutrons a não decair",
    "domando campos de 1e18 G",
    "perguntando ao ξ se está tudo bem",
    "integrando TOV com carinho",
    "desviando de buracos negros",
    "alinhando momentos anômalos",
    "chamando Tolman, Oppenheimer e Volkoff",
    "espremendo híperons",
    "medindo o raio com régua de 10 m",
    "acalmando o vácuo de Gaete-Helayël",
];

/// Barra de progresso dos estudos: um pulsar girando, a etapa, estrelas
/// enchendo a barra, tempo decorrido e restante, e uma frase que muda.
pub fn progress(total: usize, stage: &str) -> ProgressBar {
    let bar = ProgressBar::new(total as u64);
    if let Ok(style) = ProgressStyle::with_template(
        "  {spinner:.magenta} {prefix:.bold.cyan} {bar:30.yellow/blue} {pos}/{len} [{elapsed_precise}, resta ~{eta}] {msg:.dim}",
    ) {
        bar.set_style(
            style
                .tick_strings(&["✶", "✷", "✸", "✹", "✺", "✹", "✸", "✷", "★"])
                .progress_chars("★✦·"),
        );
    }
    bar.set_prefix(stage.to_string());
    bar.set_message(PROGRESS_PHRASES[0]);
    bar.enable_steady_tick(Duration::from_millis(120));
    bar
}

/// Conta uma tarefa concluída na barra e troca a frase a cada 5 s.
pub fn tick(bar: &ProgressBar) {
    bar.inc(1);
    let i = (bar.elapsed().as_secs() / 5) as usize % PROGRESS_PHRASES.len();
    bar.set_message(PROGRESS_PHRASES[i]);
}
