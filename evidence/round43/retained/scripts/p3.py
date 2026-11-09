p = "/home/ubuntu/projects/DnaCalc/DnaCalc/crates/calc-func/src/ported/regression_forecast_family.rs"
s = open(p).read()

old_doc = """/// 1. With a constant, each predictor column and `y` are centred on their means (each mean the
///    forward sum divided by the count), and a column of ones leads the matrix; without one, the
///    raw columns are the matrix.
/// 2. Householder QR. Before each predictor column's step, the predictor column with the largest
///    running remainder (5 below) comes next; on an exact tie the later column wins with a
///    constant and the earlier one without (R42). The ones column always goes first. At step"""
new_doc = """/// 1. With a constant, each predictor column and `y` are centred on their means (each mean the
///    forward sum divided by the count), and a column of ones is stored after the predictors and
///    swapped to the front, so the first predictor moves to the last place (R43); without one,
///    the raw columns are the matrix.
/// 2. Householder QR. Before each predictor column's step, the predictor column with the largest
///    running remainder (5 below, a negative one read as zero) is swapped in; on an exact tie the
///    earlier column in the current order wins (R43: with a constant this is the later one in
///    the given order, which R42 read as a constant-dependent tie). At step"""
assert old_doc in s
s = s.replace(old_doc, new_doc)

old5 = """///    is below 2^-41 times its sum of squares before the steps, or its entries from the diagonal
///    down are all zero. The decomposition then stops: every column not yet stepped is left out
///    with it (R42)."""
new5 = """///    is below 2^-41 times its sum of squares before the steps, or its entries from the diagonal
///    down are all zero. While other columns are still to come the test reads the reflector's
///    `vv` in place of that remainder (R43). The decomposition then stops: every column not yet
///    stepped is left out with it (R42)."""
assert old5 in s
s = s.replace(old5, new5)

old_r42 = """/// 1,303 for each (the previous program 903 and 932). Some zero-column layouts beside a
/// collinear pair still differ. Private QR and residual arithmetic remain incomplete."""
new_r42 = """/// 1,303 for each (the previous program 903 and 932).
/// R43 (Excel 16.0.20527, 2026-11-09): the ones column's place, the clipped tie-to-the-earlier
/// pivot and the `vv` test, found on the zero-column layouts R42 left open and on threshold maps
/// (n = 4 to 32: the switch is 2^-41 of the remainder for the last column and 2^-41 of `vv`
/// when a zero or tiny column is still to come). A fresh set frozen before asking (2,439
/// selections per function, 2,414 answered) matched 2,413 and 2,414 of the answered (the
/// previous program 2,225 and 2,255). Private QR and residual arithmetic remain incomplete."""
assert old_r42 in s
s = s.replace(old_r42, new_r42)

anchor = """    /// dc-j9k.14: with the larger column of a near-collinear pair kept and the other left out,"""
tests = """    /// dc-j9k.14 (R43): with a constant the ones column is stored after the predictors and
    /// swapped to the front, and a pivot tie goes to the earlier column of that order. For
    /// known_x {0, x, -x} the pivot is x, -x is collinear and the zero column is left out with
    /// it. Excel 16.0.20527, 2026-11-09 (R42 pivot set 3, DnaCalc-Research evidence/round42).
    #[cfg_attr(not(target_arch = "wasm32"), test)]
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test::wasm_bindgen_test)]
    fn ones_column_swapped_in_decides_a_mirrored_tie() {
        let x = [21.0, -23.0, 2.0, -24.0, -5.0, 9.0, 3.0, 15.0];
        let rows: Vec<Vec<f64>> = x.iter().map(|v| vec![0.0, *v, -*v]).collect();
        let rows: Vec<&[f64]> = rows.iter().map(|r| r.as_slice()).collect();
        let y = col(&[19.0, 6.0, 7.0, 12.0, -6.0, 3.0, 6.0, 7.0]);
        let cells = block(
            eval_linest_surface(&[y, matrix(&rows), bool_arg(true), bool_arg(true)]).unwrap(),
        );
        assert_eq!(cells[0], Ok(0)); // -x, collinear
        assert_eq!(cells[1], Ok(0x3fb874848140d471)); // x
        assert_eq!(cells[2], Ok(0)); // the zero column, left out with -x
        assert_eq!(cells[3 * 4 + 1], Ok(6.0f64.to_bits())); // degrees of freedom
        assert_eq!(cells[4 * 4 + 1], Ok(0x4075241d892ddca4)); // residual sum of squares
        assert_eq!(cells[2 * 4 + 1], Ok(0x401e089c4d08d0b2)); // standard error of y
    }

    /// dc-j9k.14 (R43): while a column is still to come the collinearity test reads the
    /// reflector's vv. Here the second column is the third halved with 2.5 * 2^-18 added to its last
    /// entry: its remainder is below 2^-41 of its sum of squares but vv is not, so all three
    /// predictors are kept over zero degrees of freedom. Excel 16.0.20527, 2026-11-09, a fresh
    /// question frozen before asking (R43 proof set, DnaCalc-Research evidence/round43).
    #[cfg_attr(not(target_arch = "wasm32"), test)]
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen_test::wasm_bindgen_test)]
    fn collinearity_test_reads_vv_while_columns_remain() {
        let rows: [[u64; 3]; 4] = [
            [0x0000000000000000, 0x4029000000000000, 0x4039000000000000],
            [0x0000000000000000, 0xc016000000000000, 0xc026000000000000],
            [0x0000000000000000, 0x4022000000000000, 0x4032000000000000],
            [0x0000000000000000, 0x4014000280000000, 0x4024000000000000],
        ];
        let x = CalcValue::array(
            CalcArray::from_rows(
                rows.iter()
                    .map(|r| r.iter().map(|b| CalcValue::number(f64::from_bits(*b))).collect())
                    .collect(),
            )
            .unwrap(),
        );
        let y = col(&[7.0, -14.0, 12.0, 11.0]);
        let cells = block(eval_linest_surface(&[y, x, bool_arg(true), bool_arg(true)]).unwrap());
        assert_eq!(cells[0], Ok(0xc11f4e5bd4c16924));
        assert_eq!(cells[1], Ok(0x412f4e5e87b2a7ab));
        assert_eq!(cells[2], Ok(0)); // the zero column, kept
        assert_eq!(cells[3], Ok(0xc0161f629af00000)); // intercept
        assert_eq!(cells[3 * 4 + 1], Ok(0)); // degrees of freedom
    }

"""
assert anchor in s
s = s.replace(anchor, tests + anchor)
open(p, "w").write(s)
print("ok")
