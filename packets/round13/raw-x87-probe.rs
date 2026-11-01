//! Raw x87 leaf replay for the Windows CPU comparison, from the public Intel ISA.
//! Original implementation, MIT. Compile with rustc; no product or library links.
#[cfg(not(target_arch = "x86_64"))]
compile_error!("This probe requires native x86_64 x87 instructions");

use std::{env, fs, io::Write, path::Path};

fn decode(text: &str) -> Result<[u8; 10], String> {
    if text.len() != 20 || !text.is_ascii() {
        return Err("expected twenty hex digits for ten little-endian bytes".into());
    }
    let mut result = [0; 10];
    for (i, byte) in result.iter_mut().enumerate() {
        *byte = u8::from_str_radix(&text[2 * i..2 * i + 2], 16).map_err(|e| e.to_string())?;
    }
    Ok(result)
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|x| format!("{x:02x}")).collect()
}

fn instruction(kind: &str, a: &[u8; 10], b: &[u8; 10], cw: u16) -> Result<[u8; 10], String> {
    if ![0x123f, 0x133f].contains(&cw) {
        return Err("only RN PC53/PC64 words are admitted".into());
    }
    let mut result = [0; 10];
    let mut saved = 0_u16;
    // SAFETY: the control word is restored below, all pointers name live locals,
    // and each branch consumes all of the x87 stack entries it pushes.
    unsafe {
        core::arch::asm!("fnstcw word ptr [{p}]", "fldcw word ptr [{q}]",
            p=in(reg)&mut saved, q=in(reg)&cw, options(nostack, preserves_flags));
        match kind {
            "fyl2x" => core::arch::asm!(
                "fld tbyte ptr [{a}]", "fld tbyte ptr [{b}]", "fyl2x", "fstp tbyte ptr [{o}]",
                a=in(reg)a.as_ptr(), b=in(reg)b.as_ptr(), o=in(reg)result.as_mut_ptr(),
                options(nostack,preserves_flags)),
            "fyl2xp1" => core::arch::asm!(
                "fld tbyte ptr [{a}]", "fld tbyte ptr [{b}]", "fyl2xp1", "fstp tbyte ptr [{o}]",
                a=in(reg)a.as_ptr(), b=in(reg)b.as_ptr(), o=in(reg)result.as_mut_ptr(),
                options(nostack,preserves_flags)),
            "f2xm1" => core::arch::asm!(
                "fld tbyte ptr [{a}]", "f2xm1", "fstp tbyte ptr [{o}]",
                a=in(reg)a.as_ptr(), o=in(reg)result.as_mut_ptr(), options(nostack,preserves_flags)),
            "fldln2" => core::arch::asm!("fldln2", "fstp tbyte ptr [{o}]",
                o=in(reg)result.as_mut_ptr(), options(nostack,preserves_flags)),
            "fldl2e" => core::arch::asm!("fldl2e", "fstp tbyte ptr [{o}]",
                o=in(reg)result.as_mut_ptr(), options(nostack,preserves_flags)),
            _ => {
                core::arch::asm!("fldcw word ptr [{p}]", p=in(reg)&saved, options(nostack,preserves_flags));
                return Err(format!("unknown instruction {kind}"));
            }
        }
        core::arch::asm!("fldcw word ptr [{p}]", p=in(reg)&saved, options(nostack,preserves_flags));
    }
    Ok(result)
}

fn run(input: &Path, output: &Path) -> Result<(), String> {
    let requests = fs::read_to_string(input).map_err(|e| e.to_string())?;
    let mut lines = requests.lines();
    if lines.next() != Some("id\tkind\tcw\ta\tb") {
        return Err("unexpected input header".into());
    }
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)
        .map_err(|e| e.to_string())?;
    writeln!(file, "id\tkind\tcw\ta\tb\tresult").map_err(|e| e.to_string())?;
    let mut count = 0;
    for line in lines {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 5 {
            return Err("malformed input row".into());
        }
        let cw = u16::from_str_radix(fields[2], 16).map_err(|e| e.to_string())?;
        let a = decode(fields[3])?;
        let b = decode(fields[4])?;
        let result = instruction(fields[1], &a, &b, cw)?;
        writeln!(file, "{line}\t{}", hex(&result)).map_err(|e| e.to_string())?;
        count += 1;
    }
    file.sync_all().map_err(|e| e.to_string())?;
    let mut brand = Vec::new();
    if core::arch::x86_64::__cpuid(0x80000000).eax >= 0x80000004 {
        for leaf in 0x80000002..=0x80000004 {
            let c = core::arch::x86_64::__cpuid(leaf);
            for r in [c.eax, c.ebx, c.ecx, c.edx] {
                brand.extend(r.to_le_bytes());
            }
        }
    }
    println!(
        "native CPU: {}; OS: {}; arch: {}; completed leaves: {count}",
        String::from_utf8_lossy(&brand).trim_matches(char::from(0)),
        env::consts::OS,
        env::consts::ARCH
    );
    Ok(())
}

fn main() {
    let args: Vec<_> = env::args_os().collect();
    if args.len() != 3 {
        eprintln!("usage: raw-x87-probe INPUT.tsv NEW-OUTPUT.tsv");
        std::process::exit(2);
    }
    if let Err(error) = run(Path::new(&args[1]), Path::new(&args[2])) {
        eprintln!("native probe failed: {error}");
        std::process::exit(1);
    }
}
