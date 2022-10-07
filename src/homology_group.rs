use nalgebra::{DMatrix, DVector};


pub struct SimplicalComplex {
    verts: Vec<usize>,
    edges: Vec<[usize;2]>,
    faces: Vec<[usize;3]>,
}

fn perm_sgn(){

}

fn simplex_boundary(s: &[usize], c1: &Vec<[usize]>){
    let chain = DVector::<isize>::zeros(c1.len());
    for i in 0..s.len() {
        
    }
}

// dim c1 < dim c2 
fn calc_ith_boundary(c1:&Vec<[usize]>, c2:&Vec<[usize]>){
    let mut d = DMatrix::<isize>::zeros(c2.len(), c1.len());
    for j in 0..c2.len() {
        let s = c2[j];
        d[j] = simplex_boundary(s, c1);
    }
    d.transpose()
}

fn calc_boundary_operators(sc: SimplicalComplex) -> {
    let size = 3;
    let mut boundary_list = vec![[0, 0, 0]];
    for i in 0 ..size - 1 {
        boundary_list.append(calc_ith_boundary)
    }

}

fn col_r(){

}

fn col_r_same_tor(){

}

fn calc_cohomology(){

}

fn calc_ith_homology(){

}

pub fn calc_homology_group_list(sc: SimplicalComplex) {
    let boundaries = calc_boundary_operators(sc);

}