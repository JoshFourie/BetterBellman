use group::Wnaf;
use pairing::Engine;

use crate::{domain, arith};
use domain::Domain;
use arith::Group;

use super::{KeyPairAssembly, ParameterAssembly};

pub struct Windows<E>
where
    E: Engine
{
    g1: Wnaf<(), Vec<E::G1>, Vec<i64>>,
    g2: Wnaf<(), Vec<E::G2>, Vec<i64>>
}

impl<E> Default for Windows<E>
where
    E: Engine
{
    fn default() -> Self {
        Windows {
            g1: Wnaf::new(),
            g2: Wnaf::new()
        }
    }
}

impl<E> Windows<E> 
where
    E: Engine
{
    pub fn as_based<'a,C,G>(&'a mut self, kp: &KeyPairAssembly<E>, param: &ParameterAssembly<E,C>, dom: &Domain<E,G>) -> BasedWindows<'a,E> 
    where
        G: Group<'a,E>
    {
        let domain_size: usize = dom.as_ref().len() - 1;
        BasedWindows::new(self, kp, param, domain_size)
    }
}

pub struct BasedWindows<'a,E>
where
    E: Engine
{
    pub g1: Wnaf<usize, &'a [E::G1], &'a mut Vec<i64>>,
    pub g2: Wnaf<usize, &'a [E::G2], &'a mut Vec<i64>>
}

impl<'a,E> BasedWindows<'a,E>
where
    E: Engine,
{
    fn new<C>(wind: &'a mut Windows<E>, kp: &KeyPairAssembly<E>, param: &ParameterAssembly<E,C>, domain_size: usize) -> Self {

        let (g1_query, g2_query): _ = get_queries(kp, domain_size);

        let based_g1: Wnaf<_, &'a _, &'a mut _> = wind.g1.base(param.groups.g1, g1_query);
        let based_g2: Wnaf<_, &'a _, &'a mut _> = wind.g2.base(param.groups.g2, g2_query);
            
        BasedWindows {
            g1: based_g1,
            g2: based_g2
        }
    }
}

fn get_queries<E>(kp: &KeyPairAssembly<E>, domain_size: usize) -> (usize, usize) 
where
    E: Engine
{
    let g1_query: usize = {
        let h_query: _ = domain_size;
        let icl_query: _ = kp.num.inputs + kp.num.aux;
        let a_query: _ = kp.num.inputs + kp.num.aux;
        let b_query: _ = kp.num.inputs + kp.num.aux;
        h_query + icl_query + a_query + b_query
    };
    let g2_query: usize = kp.num.inputs + kp.num.aux;

    (g1_query, g2_query)
}
