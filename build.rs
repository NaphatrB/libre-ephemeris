// MIT License
//
// Copyright (c) 2026 Libre Ephemeris Contributors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

fn main() {
    let crate_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let config = cbindgen::Config {
        language: cbindgen::Language::C,
        ..Default::default()
    };
    cbindgen::Builder::new()
        .with_crate(crate_dir.clone())
        .with_config(config)
        .generate()
        .expect("cbindgen failed")
        .write_to_file("include/le_ephemeris.h");

    generate_nutation_tables();
}

fn generate_nutation_tables() {
    let tab5_3a = std::path::Path::new("data").join("tab5.3a.txt");
    let tab5_3b = std::path::Path::new("data").join("tab5.3b.txt");
    let out = std::path::Path::new(&std::env::var("OUT_DIR").unwrap()).join("nutation_terms.rs");

    if !tab5_3a.exists() || !tab5_3b.exists() {
        println!("cargo:warning=IERS tables not found; using 77-term subset");
        write_placeholder(&out);
        return;
    }
    let a = std::fs::read_to_string(&tab5_3a).unwrap();
    let b = std::fs::read_to_string(&tab5_3b).unwrap();

    let (terms_a_j0, terms_a_j1) = parse_block_dual(&a);
    let (terms_b_j0, terms_b_j1) = parse_block_dual(&b);

    if terms_a_j0.is_empty() || terms_b_j0.is_empty() {
        println!("cargo:warning=Failed to parse IERS tables (j0: {} + {} terms)", terms_a_j0.len(), terms_b_j0.len());
        write_placeholder(&out);
        return;
    }

    // Merge j=0 terms by line index
    let n0 = terms_a_j0.len().min(terms_b_j0.len());
    // Merge j=1 terms by line index
    let n1 = terms_a_j1.len().min(terms_b_j1.len());

    let mut code = String::from(
        "// Generated from IERS Conventions 2010 Tables 5.3a/5.3b (IAU 2000_R06)\n"
    );
    code.push_str("// Unit: 0.1 microarcsecond\n");
    code.push_str("#[allow(clippy::unreadable_literal)]\n");
    code.push_str("pub const NUT_TERMS_FULL: &[(i32,i32,i32,i32,i32,i32,i32,i32,i32,i32,i32,i32,i32,i32,f64,f64,f64,f64,f64,f64,f64,f64)] = &[\n");

    // j=0 terms (constant coefficients)
    for i in 0..n0 {
        let ta = &terms_a_j0[i];
        let tb = &terms_b_j0[i];
        // tab5.3a: A_i (sine), A"_i (cosine)
        // tab5.3b: B"_i (sine), B_i (cosine) — note swapped order
        let a_i = ta.coeffs[0] * 10.0;
        let a2_i = ta.coeffs[1] * 10.0;
        let b2_i = tb.coeffs[0] * 10.0; // B"_i
        let b_i = tb.coeffs[1] * 10.0;  // B_i
        code.push_str(&format!(
            "({},{},{},{},{},{},{},{},{},{},{},{},{},{},{:.1},{:.1},{:.1},{:.1},0.0,0.0,0.0,0.0),\n",
            ta.mults[0], ta.mults[1], ta.mults[2], ta.mults[3], ta.mults[4],
            ta.mults[5], ta.mults[6], ta.mults[7], ta.mults[8], ta.mults[9],
            ta.mults[10], ta.mults[11], ta.mults[12], ta.mults[13],
            a_i, a2_i, b_i, b2_i,
        ));
    }

    // j=1 terms (time-dependent coefficients)
    for i in 0..n1 {
        let ta = &terms_a_j1[i];
        let tb = &terms_b_j1[i];
        let a_i = 0.0;
        let a2_i = 0.0;
        let b_i = 0.0;
        let b2_i = 0.0;
        // tab5.3a: A'_i (sine * t), A"'_i (cosine * t)
        // tab5.3b: B"'_i (sine * t), B'_i (cosine * t) — swapped order
        let at_i = ta.coeffs[0] * 10.0;
        let at2_i = ta.coeffs[1] * 10.0;
        let bt2_i = tb.coeffs[0] * 10.0;
        let bt_i = tb.coeffs[1] * 10.0;
        code.push_str(&format!(
            "({},{},{},{},{},{},{},{},{},{},{},{},{},{},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1},{:.1}),\n",
            ta.mults[0], ta.mults[1], ta.mults[2], ta.mults[3], ta.mults[4],
            ta.mults[5], ta.mults[6], ta.mults[7], ta.mults[8], ta.mults[9],
            ta.mults[10], ta.mults[11], ta.mults[12], ta.mults[13],
            a_i, a2_i, b_i, b2_i,
            at_i, at2_i, bt_i, bt2_i,
        ));
    }

    code.push_str("];\n");

    std::fs::write(&out, code).unwrap();
    println!("cargo:warning=Generated nutation: {} j0 + {} j1 = {} full terms", n0, n1, n0 + n1);
}

struct TermData {
    mults: Vec<i32>,
    coeffs: Vec<f64>,
}

/// Parse all data blocks from an IERS table file.
/// Returns (j0_terms, j1_terms).
fn parse_block_dual(content: &str) -> (Vec<TermData>, Vec<TermData>) {
    let lines: Vec<&str> = content.lines().collect();
    let mut j0 = Vec::new();
    let mut j1 = Vec::new();
    let mut section: i32 = 0;

    let mut in_data = false;

    for line in &lines {
        let trimmed = line.trim();

        if trimmed.is_empty() {
            continue;
        }

        // Detect section header: "j = 0" or "j = 1"
        if trimmed.starts_with("j =") || trimmed.starts_with("j =") {
            if let Some(pos) = trimmed.find('=') {
                let rest = trimmed[pos+1..].trim();
                if let Some(first) = rest.split_whitespace().next() {
                    if let Ok(s) = first.parse::<i32>() {
                        section = s;
                        in_data = false;
                        continue;
                    }
                }
            }
        }

        // Skip table headers and separators
        if trimmed.starts_with("---") || trimmed.starts_with("Table") 
            || trimmed.starts_with('i') || trimmed.starts_with("(unit")
            || trimmed.starts_with("Sum_") || trimmed.starts_with("The T")
            || trimmed.starts_with("The e") || trimmed.starts_with("The f")
        {
            in_data = false;
            continue;
        }

        // Detect data line: starts with a positive integer
        if !in_data {
            if let Some(first) = trimmed.split_whitespace().next() {
                if let Ok(idx) = first.parse::<usize>() {
                    if idx > 0 && idx <= 2000 {
                        in_data = true;
                    } else {
                        continue;
                    }
                } else {
                    continue;
                }
            } else {
                continue;
            }
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 17 { continue; }

        let idx: usize = parts[0].parse().unwrap_or(0);
        if idx == 0 { continue; }

        let mut coeffs = Vec::new();
        let mut mults = Vec::new();

        if let Ok(v) = parts[1].parse::<f64>() { coeffs.push(v); }
        else { continue; }
        if let Ok(v) = parts[2].parse::<f64>() { coeffs.push(v); }
        else { continue; }

        for p in parts.iter().skip(3).take(14) {
            if let Ok(v) = p.parse::<i32>() { mults.push(v); }
            else { break; }
        }
        if mults.len() < 14 { continue; }

        match section {
            0 => j0.push(TermData { mults, coeffs }),
            1 => j1.push(TermData { mults, coeffs }),
            _ => {}
        }
    }

    (j0, j1)
}

fn write_placeholder(path: &std::path::Path) {
    std::fs::write(path, b"pub const NUT_TERMS_FULL: &[(i32,i32,i32,i32,i32,i32,i32,i32,i32,i32,i32,i32,i32,i32,f64,f64,f64,f64,f64,f64,f64,f64)] = &[];").unwrap();
}
