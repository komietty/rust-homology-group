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

fn calc_boundary_operators(sc: SimplicalComplex) -> (DMatrix<isize>, DMatrix<isize>, DMatrix<isize>, DMatrix<isize>){
    let d0 = DMatrix::<isize>::zeros(1, sc.verts.len());
    let d1 = calc_ith_boundary(&sc.verts, &sc.edges);
    let d2 = calc_ith_boundary(&sc.edges, &sc.faces);
    let d3 = DMatrix::<isize>::from_element(0, 0, 0);
    (d0, d1, d2, d3)
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
    let b = u.slice((0, 0), (u.nrows(), r)).clone_owned();
    let t = o.b.slice((0, 0), (r, r)).clone().diagonal().clone();
    //println!("U: {}", u);
    //println!("B: {}", b);
    //println!("Z: {}", z);
    //println!("--------");
    //println!("T: {}", t);
    let zz = col_r(z.into_owned());
    let (bb, tt) = col_r_same_tor(cast_f2i(b), t);
    (bb, tt, zz)
}

fn calc_ith_homology(d1: DMatrix<isize>, d2: DMatrix<isize>) -> (DMatrix<i128>, Vec<i128>){
    //println!("d1: {}", d1.clone());
    let mut z1 = DMatrix::<i128>::identity(d1.ncols(), d1.ncols());
    let mut b0 = DMatrix::<i128>::zeros(1, 1);
    let mut t0 = DMatrix::<i128>::zeros(1, 1);
    let mut b1 = DMatrix::<i128>::zeros(1, 1);
    let mut t1 = TMatrix::from_row_slice(&[]);
    if !is_zeros(&d1) {
        let set = calc_cohomology(d1);
        b0 = set.0;
        z1 = set.2;
    }

    if d2.is_empty() {
        b1 = DMatrix::<i128>::from_element(0, 0, 0);
    } else {
        let set = calc_cohomology(d2);
        b1 = set.0;
        t1 = set.1;
    }
    let z_nonzero = np_map(&z1.transpose());
    let b_nonzero = np_map(&b1.transpose()); 
    let mut non_trivial = vec![];
    //for v in z_nonzero { println!("znonzero: {}", v); }
    //for v in b_nonzero { println!("bnonzero: {}", v); }

    for (i, z) in z_nonzero.clone().iter().enumerate() { 
        if !b_nonzero.contains(z) {
            non_trivial.push((i, 1));
        }
    }

    for (i, b) in b_nonzero.iter().enumerate() { 
        if t1[(i, 0)] > 1 {
            non_trivial.push((z_nonzero.clone().into_iter().find(|z| z == b).unwrap(), t1[i]));
        }
    }
    //println!("z1: {}", z1);
    //println!("b1: {}", b1);

    //print!("nt:");
    //for k in non_trivial.clone() { print!("({}, {})", k.0, k.1); }
    //println!("");

    let mut zz = DMatrix::<i128>::zeros( z1.nrows(), non_trivial.len());
    let mut tt = vec![];
    for (i, x) in non_trivial.iter().enumerate() {
        zz.set_column(i, &z1.column(x.0));
        tt.push(x.1);
    };

    println!("zz: {}", zz);
    //print!("z1: {}", z1);
    print!("k: ");
    for k in tt.clone() { print!("{}, ", k); }
    println!();
    println!("-----");

    (zz, tt)
}

pub fn calc_homology_group_list(sc: SimplicalComplex) {
    let (d0, d1, d2, d3) = calc_boundary_operators(sc);
    let mut homology_groups = vec![];
    homology_groups.push(calc_ith_homology(d0.clone(), d1.clone()));
    homology_groups.push(calc_ith_homology(d1.clone(), d2.clone()));
    homology_groups.push(calc_ith_homology(d2.clone(), d3.clone()));

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

#[test]
fn test_2() {
    let sc = SimplicalComplex{
        verts: vec![vec![0], vec![1], vec![2], vec![3]],
        edges: vec![
            vec![0, 1], vec![0, 2], vec![0, 3],
            vec![1, 2], vec![1, 3], vec![2, 3],
            ],
        faces: vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3]],
    };
    calc_homology_group_list(sc);
}

#[test]
fn test_3() {
    let sc = SimplicalComplex{
        verts: vec![
            vec![0],
            vec![1],
            vec![2],
            vec![3],
            vec![4],
            vec![5],
            vec![6],
            vec![7],
            vec![8],
            ],
        edges: vec![
            vec![0, 1],
            vec![0, 2],
            vec![0, 3],
            vec![0, 4],
            vec![0, 6],
            vec![0, 7],
            vec![1, 2],
            vec![1, 4],
            vec![1, 5],
            vec![1, 7],
            vec![2, 3],
            vec![2, 5],
            vec![2, 6],
            vec![2, 8],
            vec![3, 4],
            vec![3, 5],
            vec![3, 7],
            vec![3, 8],
            vec![4, 5],
            vec![4, 6],
            vec![4, 8],
            vec![5, 6],
            vec![5, 7],
            vec![6, 7],
            vec![6, 8],
            vec![7, 8],
            ],
        faces: vec![
            vec![0, 1, 4],
            vec![1, 4, 5],
            vec![1, 2, 5],
            vec![2, 5, 6],
            vec![0, 6, 2],
            vec![0, 6, 4],
            vec![3, 4, 5],
            vec![3, 7, 5],
            vec![5, 6, 7],
            vec![6, 7, 8],
            vec![4, 8, 6],
            vec![3, 4, 8],
            vec![0, 3, 7],
            vec![0, 1, 7],
            vec![1, 7, 8],
            vec![1, 2, 8],
            vec![2, 8, 3],
            vec![0, 3, 2]
        ],
    };
    calc_homology_group_list(sc);
}

#[test]
fn test_4() {
    let sc = SimplicalComplex{
        verts: vec![
            vec![0],
            vec![1],
            vec![2],
            vec![3],
            vec![4],
            vec![5],
            vec![6],
            vec![7],
            vec![8],
            ],
        edges: vec![
            vec![0, 1],
            vec![0, 2],
            vec![0, 3],
            vec![0, 4],
            vec![0, 6],
            vec![0, 7],
            vec![1, 2],
            vec![1, 5],
            vec![1, 6],
            vec![1, 7],
            vec![2, 3],
            vec![2, 4],
            vec![2, 5],
            vec![2, 8],
            vec![3, 4],
            vec![3, 5],
            vec![3, 7],
            vec![3, 8],
            vec![4, 5],
            vec![4, 6],
            vec![4, 8],
            vec![5, 6],
            vec![5, 7],
            vec![6, 7],
            vec![6, 8],
            vec![7, 8],
            ],
        faces: vec![
            vec![0, 4, 2],
            vec![2, 4, 5],
            vec![1, 2, 5],
            vec![1, 5, 6],
            vec![0, 1, 6],
            vec![0, 6, 4],
            vec![3, 5, 4],
            vec![3, 7, 5],
            vec![5, 7, 6],
            vec![6, 7, 8],
            vec![4, 6, 8],
            vec![3, 4, 8],
            vec![0, 7, 3],
            vec![0, 1, 7],
            vec![1, 8, 7],
            vec![1, 2, 8],
            vec![2, 3, 8],
            vec![0, 3, 2]
        ],
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