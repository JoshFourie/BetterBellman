use pairing::Engine;

pub struct KeyPairNum {
    pub inputs: usize,
    pub aux: usize,
    pub constraints: usize
}

impl Default for KeyPairNum {
    fn default() -> Self {
        Self {
            inputs: 0,
            aux: 0,
            constraints: 0
        }
    }
}

pub struct KeyPairWires<E> 
where
    E: Engine
{
    pub at: Vec<Vec<(E::Fr, usize)>>,
    pub bt: Vec<Vec<(E::Fr, usize)>>,
    pub ct: Vec<Vec<(E::Fr, usize)>>,
}

impl<E> KeyPairWires<E>
where
    E: Engine
{
    pub fn flatten(&self) -> FlatKeyPairWires<'_,E> {
        FlatKeyPairWires::from(self)
    }
}


impl<E> Default for KeyPairWires<E>
where
    E: Engine
{
    fn default() -> Self {
        Self {
            at: Vec::new(), 
            bt: Vec::new(),
            ct: Vec::new(),
        }
    }
}

pub struct FlatKeyPairWires<'a, E: Engine>(Vec<(&'a [(E::Fr, usize)], &'a [(E::Fr, usize)], &'a [(E::Fr, usize)])>);

impl<'a,E> FlatKeyPairWires<'a,E> 
where
    E: Engine
{
    pub fn chunks(&self, chunk_size: usize) -> std::slice::Chunks<'_, (&[(E::Fr, usize)], &[(E::Fr, usize)], &[(E::Fr, usize)])> {
        self.0.chunks(chunk_size)
    }
}

impl<'a,E> From <&'a KeyPairWires<E>> for FlatKeyPairWires<'a,E> 
where
    E: Engine
{
    fn from(kp: &'a KeyPairWires<E>) -> Self {
        let flattened: Vec<_> = kp.at.iter()
            .zip(kp.bt.iter())
            .zip(kp.ct.iter())
            .map(|((at, bt), ct)| {
                (at.as_slice(), bt.as_slice(), ct.as_slice())
            }).collect();
        FlatKeyPairWires(flattened)
    }
}