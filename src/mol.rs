use crate::{atom::Atom, bond::BondWithMol};
use petgraph::graph::{Graph, NodeIndex};

#[derive(Debug)]
pub struct Mol {
    graph: Graph<Atom, BondWithMol>,
}

impl Mol {
    pub fn atom_with_idx(&self, idx: usize) -> Option<&Atom> {
        self.graph.node_weight(NodeIndex::new(idx))
    }

    pub fn atom_bonds(&self, atom: &Atom) -> impl Iterator<Item = &BondWithMol> {
        let node_idx = NodeIndex::new(atom.index);
        self.graph.edges(node_idx).map(|edge| edge.weight())
    }

    pub fn atom_degree(&self, atom: &Atom) -> i32 {
        let node_idx = NodeIndex::new(atom.index);
        self.graph.edges(node_idx).count() as i32 // port: TODO review i32 cast
    }

    pub fn atom_neighbors(&self, atom: &Atom) -> impl Iterator<Item = &Atom> {
        let node_idx = NodeIndex::new(atom.index);
        self.graph
            .neighbors(node_idx)
            .filter_map(move |neighbor_idx| self.graph.node_weight(neighbor_idx))
    }
}
