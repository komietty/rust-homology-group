//
// smith normal form with nalgebra
//

use nalgebra::{DMatrix, DMatrixSlice, OMatrix, Dynamic};
type IM = OMatrix<i128, Dynamic, Dynamic>;

pub struct Decomposed {
    pub p: IM,
    pub q: IM,
    pub b: IM,
}

fn iamin_full(m: DMatrixSlice<i128>) -> (usize, usize) {
    assert!(!m.is_empty(), "The input matrix must not be empty.");
    let mut min = unsafe { m.get_unchecked((0, 0)).abs() };
    let mut idx = (0, 0);
    for j in 0..m.ncols() {
        for i in 0..m.nrows() {
            let val = unsafe { m.get_unchecked((i, j)).abs() };
            if (val < min || min == 0) && val > 0 {
                min = val;
                idx = (i, j);
            }
        }
    }
    idx
}

fn is_mod_zeros(m: DMatrixSlice<i128>, v: i128) -> bool {
    let (nr, nc) = m.shape();
    for j in 0..nc {
        for i in 0..nr {
            let e = unsafe { m.get_unchecked((i, j)) };
            if e % v != 0 && e.abs() > v {
                return false;
            }
        }
    }
    true
}

fn eij(i: usize, j: usize, n: usize) -> IM {
    let mut m = DMatrix::identity(n, n);
    m[(i, i)] = 0;
    m[(j, j)] = 0;
    m[(i, j)] = 1;
    m[(j, i)] = 1;
    m
}

fn ei(i: usize, n: usize) -> IM {
    let mut m = DMatrix::identity(n, n);
    m[(i, i)] = -1;
    m
}

fn ec(i: usize, j: usize, c: i128, n: usize) -> IM {
    let mut m = DMatrix::identity(n, n);
    m[(i, j)] = c;
    m
}

fn row_null(a: &IM, k: usize) -> Decomposed {
    let mut b = a.clone();
    let mut p = IM::identity(a.nrows(), a.nrows());
    let q = IM::identity(a.ncols(), a.ncols());
    for i in k + 1..a.nrows() {
        let d = a[(i, k)] / a[(k, k)];
        let e = ec(i, k, -d, a.nrows());
        b = &e * &b;
        p = &e * &p;
    }
    Decomposed { p: p, q: q, b: b }
}

fn col_null(a: &IM, k: usize) -> Decomposed {
    let mut b = a.clone();
    let mut q = IM::identity(a.ncols(), a.ncols());
    let p = IM::identity(a.nrows(), a.nrows());
    for j in k + 1..a.ncols() {
        let d = a[(k, j)] / a[(k, k)];
        let e = ec(k, j, -d, a.ncols());
        b = &b * &e;
        q = &q * &e;
    }
    Decomposed { p: p, q: q, b: b }
}

fn mod_full(m: DMatrixSlice<i128>, val: i128) -> (usize, usize) {
    assert!(!m.is_empty(), "The input matrix must not be empty.");
    let mut idx = (0, 0);
    'l: for j in 0..m.ncols() {
        for i in 0..m.nrows() {
            let e = unsafe { m.get_unchecked((i, j)) };
            if e % val != 0 {
                idx = (i, j);
                break 'l;
            }
        }
    }
    idx
}

fn swap_min(a: &IM, k: usize) -> Decomposed {
    let s = a.slice((k, k), (a.nrows() - k, a.ncols() - k));
    let m = iamin_full(s);
    let (i, j) = (m.0 + k, m.1 + k);
    let p = eij(k, i, a.nrows());
    let mut q = eij(j, k, a.ncols());
    let mut b = &p * a * &q;
    if b[(k, k)] < 0 {
        let i = ei(k, a.ncols());
        b *= &i;
        q *= &i;
    }
    Decomposed { p: p, q: q, b: b }
}


fn remn_mod(a: &IM, k: usize) -> Decomposed {
    let (nr, nc) = a.shape();
    let mut p = IM::identity(nr, nr);
    let mut q = IM::identity(nc, nc);
    let idx = mod_full(a.slice((k + 1, k + 1), (nr - k - 1, nc - k - 1)), a[(k, k)]);
    let (i, j) = (idx.0 + k + 1, idx.1 + k + 1);
    let d = a[(i, j)] / a[(k, k)];
    let i1 = ec(i, k, d, nr);
    let i2 = ec(k, j, -1, nc);
    let b = &i1 * a * &i2;
    p *= i1;
    q *= i2;
    Decomposed { p: p, q: q, b: b }
}

pub fn smith_normalize(a: &IM) -> Decomposed {
    let iszeros    = |m: DMatrixSlice<i128>| m.iter().all(|v| *v == 0);
    let (nr, nc) = a.shape();
    let mut b = a.clone();
    let mut p = IM::identity(nr, nr);
    let mut q = IM::identity(nc, nc);
    for k in 0..if nr < nc { nr } else { nc } {
        if iszeros(b.slice((k, k), (nr - k, nc - k))) {
            break;
        }
        loop {
            let d1 = swap_min(&b, k);
            let d2 = row_null(&d1.b, k);
            (p, q) = (&d2.p * &d1.p * &p, &q * &d1.q * &d2.q);
            if iszeros(d2.b.slice((k + 1, k), (nr - k - 1, 1))) {
                let d3 = col_null(&d2.b, k);
                (p, q) = (&d3.p * &p, &q * &d3.q);
                if iszeros(d3.b.slice((k, k + 1), (1, nc - k - 1))) {
                    if is_mod_zeros(
                        d3.b.slice((k + 1, k + 1), (nr - k - 1, nc - k - 1)),
                        d3.b[(k, k)],
                    ) {
                        b = d3.b;
                        break;
                    } else {
                        let d4 = remn_mod(&d3.b, k);
                        (p, q, b) = (&d4.p * &p, &q * &d4.q, d4.b.clone());
                    }
                } else {
                    b = d3.b;
                }
            } else {
                b = d2.b;
            }
        }
    }
    Decomposed { p: p, q: q, b: b }
}


#[test]
fn snf_test() {
    let m = DMatrix::from_row_slice( 4, 4, &[
            -6, 111, -36, 6, 5, -672, 210, 74, 0, -255, 81, 24, -7, 255, -81, -10,
        ],
    );
    assert_eq!(
        smith_normalize(&m).b,
        DMatrix::from_row_slice(4, 4, &[1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 21, 0, 0, 0, 0, 0,])
    );
}
