// instead of cor opatation with torsion, rref might be simpler 

use nalgebra::{DMatrix, DVector, DMatrixSlice, U1, Dynamic, OMatrix};
use num::integer::gcd;
use crate::snf::{smith_normalize};
pub struct SimplicalComplex {
    verts: Vec<Vec<usize>>,
    edges: Vec<Vec<usize>>,
    faces: Vec<Vec<usize>>,
}
type TM = OMatrix<i128, Dynamic, U1>;
type IM = OMatrix<i128, Dynamic, Dynamic>;
type IV = DVector<i128>;

fn perm_sign(l1: &Vec<usize>, l2: &Vec<usize>, l: usize) -> i128 {
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

fn simplex_boundary(s: &Vec<usize>, c: &Vec<Vec<usize>>) -> IV {
    let mut chain = IV::zeros(c.len());
    for i in 0..s.len() {
        let mut sb = vec![];
        let mut id = 0;
        for (j, v) in s.iter().enumerate() { if i != j { sb.push(*v); } }
        for (j, s) in c.iter().enumerate() {
            let f1 = sb.len() == s.len();
            let f2 = s.iter().all(|v| sb.contains(v));
            if f1 && f2 { id = j; }
        }
        chain[id] = (-1_i128).pow(i as u32) * perm_sign(&c[id], &sb, sb.len());
    }
    chain
}

fn calc_ith_boundary(c1:&Vec<Vec<usize>>, c2:&Vec<Vec<usize>>) -> IM {
    let mut d = IM::zeros(c2.len(), c1.len());
    for (j, s) in c2.iter().enumerate() {
        let b = simplex_boundary(&s, c1);
        d.set_row(j, &b.transpose());
    }
    d.transpose()
}

fn calc_boundary_operators(c: SimplicalComplex) -> (IM, IM, IM, IM){
    let d0 = IM::zeros(1, c.verts.len());
    let d3 = IM::from_element(0, 0, 0);
    let d1 = calc_ith_boundary(&c.verts, &c.edges);
    let d2 = calc_ith_boundary(&c.edges, &c.faces);
    (d0, d1, d2, d3)
}

fn col_r(mm: IM) -> IM {
    let mut m = mm.clone(); 
    if m.is_empty() { return m; }
    let nc = m.ncols();
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

fn col_r_same_tor(mm: IM, tt: TM) -> (IM, TM){
    let mut m = mm.clone(); 
    let mut t = tt.clone();
    if m.is_empty() { return (m, t); }
    let nc = m.ncols();
    for j1 in 0..nc {
        let i = np_map(&m.transpose())[j1];
        for j2 in 0..nc {
            let g = gcd(t[j1], t[j2]);
            let f = (t[j2] > t[j1]) || (t[j1] > 1 && t[j2] > 1 && g == 1);
            if j2 == j1 || f { continue; }
            let c = sign(m[(i, j2)]) * (m[(i, j2)].abs() / m[(i, j1)]);
            m.set_column(j2, &(m.column(j2) - m.column(j1) * c));
            t[(j2, 0)] = g;
        }
    }
    (m, t)
}

fn calc_cohomology(d: IM) -> (IM, TM, IM) {
    let r = d.clone().cast::<f64>().rank(1e-8);
    let o = smith_normalize(&d.cast::<i128>());
    let u = o.p.cast::<f64>().try_inverse().unwrap();
    let z = o.q.slice((0, r), (o.q.nrows(), o.q.ncols() - r)).clone();
    let t = o.b.slice((0, 0), (r, r)).clone().diagonal().clone();
    let b = u.slice((0, 0), (u.nrows(), r)).clone_owned();
    let zz = col_r(z.into_owned());
    let (bb, tt) = col_r_same_tor(cast_f2i(b), t);
    (bb, tt, zz)
}

fn calc_ith_homology(d1: IM, d2: IM) -> (IM, Vec<i128>){
    let mut z1 = IM::identity(d1.ncols(), d1.ncols());
    let mut b1 = IM::zeros(1, 1);
    let mut t1 = TM::from_row_slice(&[]);

    if !is_zeros(&d1) {
        let set = calc_cohomology(d1);
        z1 = set.2;
    }

    if d2.is_empty() {
        b1 = IM::from_element(0, 0, 0);
    } else {
        let set = calc_cohomology(d2);
        b1 = set.0;
        t1 = set.1;
    }

    let z_nonzero = np_map(&z1.transpose());
    let b_nonzero = np_map(&b1.transpose()); 
    let mut non_trivial = vec![];

    for (i, z) in z_nonzero.iter().enumerate() { 
        if !b_nonzero.contains(z) {
            non_trivial.push((i, 1));
        }
    }

    for (i, b) in b_nonzero.iter().enumerate() { 
        if t1[(i, 0)] > 1 {
            let mut basis = 0;
            for (j, z) in z_nonzero.iter().enumerate() { if z == b { basis = j;break; } }
            non_trivial.push((basis, t1[(i, 0)]));
        }
    }

    let mut zz = IM::zeros(z1.nrows(), non_trivial.len());
    let mut tt = vec![];
    for (i, x) in non_trivial.iter().enumerate() {
        zz.set_column(i, &z1.column(x.0));
        tt.push(x.1);
    };

    println!("-----");
    println!("zz: {}", zz);
    print!("k: ");
    for k in tt.clone() { print!("{}, ", k); }
    println!("");
    println!("-----");

    (zz, tt)
}

pub fn calc_homology_group_list(sc: SimplicalComplex) -> Vec<(IM, Vec<i128>)>{
    let (d0, d1, d2, d3) = calc_boundary_operators(sc);
    let mut g = vec![];
    g.push(calc_ith_homology(d0.clone(), d1.clone()));
    g.push(calc_ith_homology(d1.clone(), d2.clone()));
    g.push(calc_ith_homology(d2.clone(), d3.clone()));
    g
}

#[test]
fn test_2() {
    let sc = SimplicalComplex{
        verts: vec![vec![0], vec![1], vec![2], vec![3]],
        edges: vec![vec![0, 1], vec![0, 2], vec![0, 3], vec![1, 2], vec![1, 3], vec![2, 3]],
        faces: vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3]],
    };
    let rs = calc_homology_group_list(sc);
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
            vec![1, 8],
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
            vec![1, 8],
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

fn is_zeros(m: &IM) -> bool {
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

fn np_map(d: &IM) -> Vec<usize> {
    let mut v = vec![];
    for i in 0..d.nrows() {
    for j in 0..d.ncols() {
        if d[(i, j)] != 0 { v.push(j); break; }
    }
    }
    v
}

fn cast_f2i(d: DMatrix<f64>) -> IM {
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