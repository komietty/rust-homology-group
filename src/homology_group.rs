use nalgebra::{DMatrix, DVector, DMatrixSlice, U1, Dynamic, OMatrix};
use num::integer::gcd;
use crate::snf::{smith_normalize};
pub struct SimplicalComplex {
    verts: Vec<Vec<usize>>,
    edges: Vec<Vec<usize>>,
    faces: Vec<Vec<usize>>,
}

fn perm_sgn(l1: &Vec<usize>, l2: &Vec<usize>, l: usize) -> isize {
    match l {
        2 => { return if l1[0] == l2[0] { 1 } else { -1 } },
        3 => {
             let f1 = l1[0] == l2[0] && l1[1] == l2[1] && l1[2] == l2[2];
             let f2 = l1[0] == l2[1] && l1[1] == l2[2] && l1[2] == l2[0];
             let f3 = l1[0] == l2[2] && l1[1] == l2[0] && l1[2] == l2[1];
             return if f1 || f2 || f3 { 1 } else { -1 }
        }
        _ => 1
        
    }
}

fn simplex_boundary(s: &Vec<usize>, c1: &Vec<Vec<usize>>) -> DVector<isize>{
    let mut chain = DVector::<isize>::zeros(c1.len());
    for i in 0..s.len() {
        let mut sb = vec![];
        let mut id = 0;
        for (j, v) in s.iter().enumerate() { if i != j { sb.push(*v); } }
        for (j, s) in c1.iter().enumerate() {
            let f1 = sb.len() == s.len();
            let f2 = s.iter().all(|v| sb.contains(v));
            if f1 && f2 { id = j; }
        }
        chain[id] = (-1_isize).pow(i as u32) * perm_sgn(&c1[id], &sb, sb.len());
    }
    chain
}

fn calc_ith_boundary(c1:&Vec<Vec<usize>>, c2:&Vec<Vec<usize>>) -> DMatrix<isize>{
    let mut d = DMatrix::<isize>::zeros(c2.len(), c1.len());
    for (j, s) in c2.iter().enumerate() {
        let b = simplex_boundary(&s, c1);
        d.set_row(j, &b.transpose());
    }
    d.transpose()
}

fn calc_boundary_operators(sc: SimplicalComplex) -> (DMatrix<isize>, DMatrix<isize>){
    let d1 = calc_ith_boundary(&sc.verts, &sc.edges);
    let d2 = calc_ith_boundary(&sc.edges, &sc.faces);
    (d1, d2)
}

fn col_r(mm: DMatrix<i128>) -> DMatrix<i128>{
    let mut m = mm.clone(); 
    let (nr, nc) = m.shape();
    if nr * nc == 0 { return m; }
    for j1 in 0..nc {
        let i = np_map(&m.transpose())[j1];
        for j2 in 0..nc {
            if j2 == j1 { continue; }
            let c = sign(m[(i, j2)]) * (m[(i, j2)].abs() / m[(i, j1)]);
            m.set_column(j2, &(m.column(j2) - m.column(j1) * c));
        }
    }
    m
}
type TMatrix = OMatrix<i128, Dynamic, U1>;

fn col_r_same_tor(mut mm: DMatrix<i128>, tt: TMatrix)-> (DMatrix<i128>, TMatrix){
    let mut m = mm.clone(); 
    let mut t = tt.clone();
    let (nr, nc) = m.shape();
    if nr * nc == 0 { return (m, t); }
    for j1 in 0..nc {
        let i = np_map(&m.transpose())[j1];
        for j2 in 0..nc {
            let gt = gcd(t[j1], t[j2]);
            let ff = (t[j2] > t[j1]) || (t[j1] > 1 && t[j2] > 1 && gt == 1);
            if j2 == j1 || ff { continue; }
            let c = sign(m[(i, j2)]) * (m[(i, j2)].abs() / m[(i, j1)]);
            m.set_column(j2, &(m.column(j2) - m.column(j1) * c));
            t[(j2, 0)] = gt;
        }
    }
    (m, t)
}

fn calc_cohomology(d: DMatrix<isize>) -> (DMatrix<i128>, TMatrix, DMatrix<i128>) {
    let r = d.clone().cast::<f64>().rank(1e-8);
    let o = smith_normalize(&d.cast::<i128>());
    let u = o.p.cast::<f64>().try_inverse().unwrap();
    let z = o.q.slice((0, r), (o.q.nrows(), o.q.ncols() - r)).clone();
    let t = o.b.slice((0, 0), (r, r)).clone().diagonal().clone();
    let b = u.slice((0, 0), (u.nrows(), r)).clone_owned();
    //println!("U: {}", u);
    //println!("Z: {}", z);
    //println!("T: {}", t);
    //println!("B: {}", b);
    let zz = col_r(z.into_owned());
    let (bb, tt) = col_r_same_tor(cast_f2i(b), t);
    (bb, tt, zz)
}

fn calc_ith_homology(d1: DMatrix<isize>, d2: DMatrix<isize>) -> DMatrix<isize>{
    println!("d1: {}", d1.clone());
    let mut z1 = DMatrix::<i128>::identity(d1.ncols(), d1.ncols());
    let mut b0 = DMatrix::<i128>::zeros(1, 1);
    let mut t0 = DMatrix::<i128>::zeros(1, 1);
    let mut b1 = DMatrix::<i128>::zeros(1, 1);
    let mut t1 = DMatrix::<i128>::zeros(1, 1);
    if !is_zeros(&d1) {
        let set = calc_cohomology(d1);
        b0 = set.0;
        z1 = set.2;
    }
    if !d2.is_empty() {
        let set = calc_cohomology(d2);
        b1 = set.0;
    }
    let z_nonzero = np_map(&z1.transpose());
    let b_nonzero = np_map(&b1.transpose()); 
    let mut non_trivial = vec![];
    println!("z1: {}", z1);
    println!("b1: {}", b1);
    //for v in z_nonzero { println!("znonzero: {}", v); }
    //for v in b_nonzero { println!("bnonzero: {}", v); }

    for (i, z) in z_nonzero.iter().enumerate() { 
        if !b_nonzero.contains(z) {
            non_trivial.push((i, 1));
        }
    }

    for (i, b) in b_nonzero.iter().enumerate() { 
        if t1[i] > 1 {
            non_trivial.push((i, 1));
        }
    }
    DMatrix::<isize>::zeros(1, 1)
}

pub fn calc_homology_group_list(sc: SimplicalComplex) {
    let boundaries = calc_boundary_operators(sc);
    let mut homology_groups = vec![];
    let da = DMatrix::from_row_slice(1, 3, &[0,0,0]);
    let db = DMatrix::from_element(0, 0, 0);
    //print!("da:{}, b0: {}", da, boundaries.0);
    //print!("b1:{}, db: {}", boundaries.0, db);
    homology_groups.push(calc_ith_homology(da, boundaries.0.clone()));
    homology_groups.push(calc_ith_homology(boundaries.0.clone(), db));

}

#[test]
fn test() {
    let sc = SimplicalComplex{
        verts: vec![vec![0], vec![1], vec![2]],
        edges: vec![vec![0, 1], vec![1, 2], vec![2, 0]],
        faces: vec![vec![0, 1, 2]],
    };
    calc_homology_group_list(sc);
}

fn is_zeros(m: &DMatrix<isize>) -> bool {
    let (nr, nc) = m.shape();
    for j in 0..nc {
        for i in 0..nr {
            let e = unsafe { m.get_unchecked((i, j)) };
            if *e != 0 {
                return false;
            }
        }
    }
    true
}

fn np_map(d: &DMatrix<i128>) -> Vec<usize> {
    let mut v = vec![];
    for i in 0..d.nrows() {
    for j in 0..d.ncols() {
        if d[(i, j)] != 0 { v.push(j); break; }
    }
    }
    v
}

fn cast_f2i(d: DMatrix<f64>) -> DMatrix<i128> {
    let mut m = DMatrix::<i128>::zeros(d.nrows(), d.ncols());
    for i in 0..d.nrows() {
        for j in 0..d.ncols() {
            m[(i, j)] = d[(i, j)].round() as i128;
        }
    }
    m
}

fn sign(i:i128)-> i128 {
    if i == 0 { return 0; }
    else { return i / i.abs(); }
}

/*
#[test]
fn aaa() {
    let a = DMatrix::from_row_slice(3, 3, &[
        0, 0, 1,
        0, 1, 0,
        1, 0, 0,
    ]);
    let ve = np_map(&a);
    for v in ve {
        print!("{}, ", v);
    }
    println!("");
}
*/