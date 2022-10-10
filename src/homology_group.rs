use nalgebra::{DMatrix, DVector, DMatrixSlice};
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

fn col_r(m: DMatrix<i128>) -> DMatrix<i128>{
    let (nr, nc) = m.shape();
    if nr * nc == 0 { return m; }
    for j1 in 0..nc {
        let mt = m.transpose();
        let v = vec![];
        'l: for ii in 0..mt.nrows() {
            for jj in 0..mt.ncols() {
                if m[(ii, jj)] != 0 {
                    v.push(ii);
                    break 'l;
                }
            }
        }
        let i = v[j1];
        for j2 in 0..nr {
            if j2 == j1 { continue; }
            let c = m([i, j2]).sign() * (m([i, j2]) / m([i, j1]));
            m.set_column(j2, m.row(j2) - c * m[j1]);
        }
    }
    m
}

fn col_r_same_tor(mut m: DMatrix<isize>, tortion: DMatrix<isize>){

}

fn calc_cohomology(d: DMatrix<isize>) -> (DMatrix<isize>, DMatrix<isize>, DMatrix<isize>) {
    let r = d.clone().cast::<f64>().rank(1e-8);
    let o = smith_normalize(&d.cast::<i128>());
    let u = o.p.cast::<f64>().try_inverse().unwrap();
    let z = o.q.slice((0, r), (o.q.nrows(), o.q.ncols() - r)).clone();
    let t = o.b.slice((0, 0), (r, r)).clone().diagonal().clone();
    let b = u.slice((0, 0), (u.nrows(), r)).clone();
    println!("U: {}", u);
    println!("Z: {}", z);
    println!("T: {}", t);
    println!("B: {}", b);
    let z = col_r(z);
    let (b, t) = col_r_same_tor(b, t);
    (b, t, z)
}

fn calc_ith_homology(d1: DMatrix<isize>, d2: DMatrix<isize>) -> DMatrix<isize>{
    let mut z1 = DMatrix::<isize>::identity(d1.ncols(), d1.ncols());
    let mut b0 = DMatrix::<isize>::zeros(1, 1);
    let mut t0 = DMatrix::<isize>::zeros(1, 1);
    let mut b1 = DMatrix::<isize>::zeros(1, 1);
    let mut t1 = DMatrix::<isize>::zeros(1, 1);
    if !is_zeros(&d1) {
        let set = calc_cohomology(d1);
    }
    if !d2.is_empty() {
        let set = calc_cohomology(d2);
    }
    //let z_nonzero = list(map(lambda z:np.where(z!=0)[0][0],Z1.T))
    //let b_nonzero = list(map(lambda z:np.where(z!=0)[0][0],B1.T))
    DMatrix::<isize>::zeros(1, 1)
}

pub fn calc_homology_group_list(sc: SimplicalComplex) {
    let boundaries = calc_boundary_operators(sc);
    let mut homology_groups = vec![];
    let da = DMatrix::from_row_slice(3, 1, &[0,0,0]);
    let db = DMatrix::from_element(0, 0, 0);
    //print!("da:{}, b0: {}", da, boundaries.0);
    //homology_groups.push(calc_ith_homology(da, boundaries.0));
    //print!("b1:{}, db: {}", boundaries.0, db);
    homology_groups.push(calc_ith_homology(boundaries.0, db));

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

fn np_map(d: DMatrix<isize>) -> Vec<isize> {

}