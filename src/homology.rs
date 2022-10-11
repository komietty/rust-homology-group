//
// instead of cor opatation with torsion, rref might be simpler 
//
use nalgebra::{DMatrix, DVector, DMatrixSlice, U1, Dynamic, OMatrix};
use num::integer::gcd;
use crate::snf::{smith_normalize};

type TM = OMatrix<i128, Dynamic, U1>;
type IM = OMatrix<i128, Dynamic, Dynamic>;
type IV = DVector<i128>;

pub struct SimplicalComplex {
    verts: Vec<Vec<usize>>,
    edges: Vec<Vec<usize>>,
    faces: Vec<Vec<usize>>,
    //innrs: Vec<Vec<usize>>,
}

fn nonzeros(m: &IM) -> Vec<usize> {
    m.row_iter().map(|r|
        r.iter().position(|v| *v != 0).unwrap()
    ).collect()
}

fn cast_fi(d: DMatrixSlice<f64>) -> IM {
    let l: Vec<i128> = d.iter().map(|v| v.round() as i128).collect();
    DMatrix::from_column_slice(d.nrows(), d.ncols(), &l)
}

fn perm_sign(l1: &Vec<usize>, l2: &Vec<usize>) -> i128 {
    match l2.len() {
        2 => { if l1[0] == l2[0] { 1 } else { -1 } },
        3 => {
             let f1 = l1[0] == l2[0] && l1[1] == l2[1] && l1[2] == l2[2];
             let f2 = l1[0] == l2[1] && l1[1] == l2[2] && l1[2] == l2[0];
             let f3 = l1[0] == l2[2] && l1[1] == l2[0] && l1[2] == l2[1];
             if f1 || f2 || f3 { 1 } else { -1 }
        }
        _ => 1
    }
}

fn simplex_boundary(s: &Vec<usize>, c: &Vec<Vec<usize>>) -> IV {
    let mut b = IV::zeros(c.len());
    for i in 0..s.len() {
        let mut sb = vec![];
        let mut id = 0;
        for (j, v) in s.iter().enumerate() { if i != j { sb.push(*v); } }
        for (j, s) in c.iter().enumerate() {
            let f1 = sb.len() == s.len();
            let f2 = s.iter().all(|v| sb.contains(v));
            if f1 && f2 { id = j; break; }
        }
        b[id] = (-1_i128).pow(i as u32) * perm_sign(&c[id], &sb);
    }
    b
}

fn calc_boundary_operator(c1:&Vec<Vec<usize>>, c2:&Vec<Vec<usize>>) -> IM {
    let mut m = IM::zeros(c2.len(), c1.len());
    for (j, s) in c2.iter().enumerate() {
        m.set_row(j, &simplex_boundary(&s, c1).transpose());
    }
    m.transpose()
}

fn calc_boundary_operators(c: SimplicalComplex) -> (IM, IM, IM, IM){
    let m0 = IM::zeros(1, c.verts.len());
    let m3 = IM::from_element(0, 0, 0);
    let m1 = calc_boundary_operator(&c.verts, &c.edges);
    let m2 = calc_boundary_operator(&c.edges, &c.faces);
    (m0, m1, m2, m3)
}

fn col_rmv(m: &mut IM) {
    let si = |i:i128| if i == 0 { 0 } else { i / i.abs() }; 
    let nc = m.ncols();
    for j1 in 0..nc {
        let i = nonzeros(&m.transpose())[j1];
        for j2 in 0..nc {
            if j2 == j1 { continue; }
            let c = si(m[(i, j2)]) * (m[(i, j2)].abs() / m[(i, j1)]);
            m.set_column(j2, &(m.column(j2) - m.column(j1) * c));
        }
    }
}

fn col_rmv_same_tor(m: &mut IM, t: &mut TM){
    let si = |i:i128| if i == 0 { 0 } else { i / i.abs() }; 
    let nc = m.ncols();
    for j1 in 0..nc {
        let i = nonzeros(&m.transpose())[j1];
        for j2 in 0..nc {
            let g = gcd(t[j1], t[j2]);
            let f = (t[j2] > t[j1]) || (t[j1] > 1 && t[j2] > 1 && g == 1);
            if j2 == j1 || f { continue; }
            let c = si(m[(i, j2)]) * (m[(i, j2)].abs() / m[(i, j1)]);
            m.set_column(j2, &(m.column(j2) - m.column(j1) * c));
            t[j2] = g;
        }
    }
}

fn calc_cohomology(m: IM) -> (IM, TM, IM) {
    let s = smith_normalize(&m);
    let u = s.p.cast::<f64>().try_inverse().unwrap();
    let n = m.cast::<f64>().rank(1e-8);
    let (r, c) = s.q.shape();
    let mut z = s.q.slice((0, n), (r, c - n)).into_owned();
    let mut t = s.b.slice((0, 0), (n, n)).diagonal();
    let mut b = cast_fi(u.slice((0, 0), (u.nrows(), n)));
    if !z.is_empty() { col_rmv(&mut z); } 
    if !b.is_empty() { col_rmv_same_tor(&mut b, &mut t); }
    (b, t, z)
}

fn calc_homology(d1: IM, d2: IM) -> (IM, Vec<i128>){
    let iszeros = |m: &DMatrix<i128>| m.iter().all(|v| *v == 0);

    let z1 = if !iszeros(&d1) {
        calc_cohomology(d1).2
    } else {
        IM::identity(d1.ncols(), d1.ncols())
    };

    let (b1, t1) = if d2.is_empty() {
        (IM::from_element(0, 0, 0), TM::from_row_slice(&[]))
    } else {
        let o = calc_cohomology(d2);
        (o.0, o.1)
    };

    let z_nz = nonzeros(&z1.transpose());
    let b_nz = nonzeros(&b1.transpose()); 
    let mut non_trivial = vec![];

    for (i, z) in z_nz.iter().enumerate() { 
        if !b_nz.contains(z) { non_trivial.push((i, 1)); }
    }

    for (i, b) in b_nz.iter().enumerate() { 
        if t1[i] > 1 {
            let n = z_nz.iter().position(|z| z == b).unwrap();
            non_trivial.push((n, t1[i]));
        }
    }

    let mut z = IM::zeros(z1.nrows(), non_trivial.len());
    let mut t = vec![];
    for (i, x) in non_trivial.iter().enumerate() {
        z.set_column(i, &z1.column(x.0));
        t.push(x.1);
    };

    (z, t)
}

pub fn calc_homology_groups(sc: SimplicalComplex) -> Vec<(IM, Vec<i128>)>{
    let (d0, d1, d2, d3) = calc_boundary_operators(sc);
    let mut g = vec![];
    g.push(calc_homology(d0.clone(), d1.clone()));
    g.push(calc_homology(d1.clone(), d2.clone()));
    g.push(calc_homology(d2.clone(), d3.clone()));
    g
}

#[test]
fn test_torus() {
    let sc = SimplicalComplex{
        verts: vec![
            vec![0], vec![1], vec![2],
            vec![3], vec![4], vec![5],
            vec![6], vec![7], vec![8],
            ],
        edges: vec![
            vec![0, 1], vec![0, 2], vec![0, 3],
            vec![0, 4], vec![0, 6], vec![0, 7],
            vec![1, 2], vec![1, 4], vec![1, 5],
            vec![1, 7], vec![1, 8], vec![2, 3],
            vec![2, 5], vec![2, 6], vec![2, 8],
            vec![3, 4], vec![3, 5], vec![3, 7],
            vec![3, 8], vec![4, 5], vec![4, 6],
            vec![4, 8], vec![5, 6], vec![5, 7],
            vec![6, 7], vec![6, 8], vec![7, 8],
            ],
        faces: vec![
            vec![0, 1, 4], vec![1, 4, 5], vec![1, 2, 5],
            vec![2, 5, 6], vec![0, 6, 2], vec![0, 6, 4],
            vec![3, 4, 5], vec![3, 7, 5], vec![5, 6, 7],
            vec![6, 7, 8], vec![4, 8, 6], vec![3, 4, 8],
            vec![0, 3, 7], vec![0, 1, 7], vec![1, 7, 8],
            vec![1, 2, 8], vec![2, 8, 3], vec![0, 3, 2]
        ],
    };

    let res = calc_homology_groups(sc);
    let res_1st = &res[0];
    let res_2nd = &res[1];
    println!("z1_base1: {}", res_1st.0.transpose());
    print!("torsion_1: ");
    for r in res_1st.1.iter() { print!("{}, ", r); }
    println!("\n");
    println!("z2_base1: {}", res_2nd.0.column(0).transpose());
    println!("z2_base2: {}", res_2nd.0.column(1).transpose());
    print!("torsion_2: ");
    for r in res_2nd.1.iter() { print!("{}, ", r); }
    println!("\n");
}

#[test]
fn test_klein() {
    let sc = SimplicalComplex{
        verts: vec![
            vec![0], vec![1], vec![2],
            vec![3], vec![4], vec![5],
            vec![6], vec![7], vec![8],
            ],
        edges: vec![
            vec![0, 1], vec![0, 2], vec![0, 3],
            vec![0, 4], vec![0, 6], vec![0, 7],
            vec![1, 2], vec![1, 5], vec![1, 6],
            vec![1, 7], vec![1, 8], vec![2, 3],
            vec![2, 4], vec![2, 5], vec![2, 8],
            vec![3, 4], vec![3, 5], vec![3, 7],
            vec![3, 8], vec![4, 5], vec![4, 6],
            vec![4, 8], vec![5, 6], vec![5, 7],
            vec![6, 7], vec![6, 8], vec![7, 8],
            ],
        faces: vec![
            vec![0, 4, 2], vec![2, 4, 5], vec![1, 2, 5],
            vec![1, 5, 6], vec![0, 1, 6], vec![0, 6, 4],
            vec![3, 5, 4], vec![3, 7, 5], vec![5, 7, 6],
            vec![6, 7, 8], vec![4, 6, 8], vec![3, 4, 8],
            vec![0, 7, 3], vec![0, 1, 7], vec![1, 8, 7],
            vec![1, 2, 8], vec![2, 3, 8], vec![0, 3, 2]
        ],
    };

    let res = calc_homology_groups(sc);
    let res_1st = &res[0];
    let res_2nd = &res[1];
    println!("z1_base1: {}", res_1st.0.transpose());
    print!("torsion_1: ");
    for r in res_1st.1.iter() { print!("{}, ", r); }
    println!("\n");
    println!("z2_base1: {}", res_2nd.0.column(0).transpose());
    println!("z2_base2: {}", res_2nd.0.column(1).transpose());
    print!("torsion_2: ");
    for r in res_2nd.1.iter() { print!("{}, ", r); }
    println!("\n");
}
