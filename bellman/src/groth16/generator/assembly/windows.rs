use group::Wnaf;
use pairing::Engine;

use crate::{domain, error, Circuit};
use domain::{Domain, Group};
use error::Result;

use super::Assembly;

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
    pub fn as_based<'a,C,G>(&'a mut self, assembly: &Assembly<E,C>, domain: &Domain<E,G>) -> Result<BasedWindows<'a,E>>
    where
        G: Group<'a,E>,
        C: Circuit<E>
    {
        let domain_size: usize = domain.as_ref().len() - 1;
        BasedWindows::new(self, assembly, domain_size)
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
    fn new<C>(wind: &'a mut Windows<E>, assembly: &Assembly<E,C>, domain_size: usize) -> Result<Self> 
    where
        C: Circuit<E>
    {

        let (g1_query, g2_query): _ = get_queries(&assembly, domain_size)?;

        let based_g1: Wnaf<_, &'a _, &'a mut _> = wind.g1.base(assembly.param.as_ref()?.groups.g1, g1_query);
        let based_g2: Wnaf<_, &'a _, &'a mut _> = wind.g2.base(assembly.param.as_ref()?.groups.g2, g2_query);
            
        Ok(BasedWindows {
            g1: based_g1,
            g2: based_g2
        })
    }
}

fn get_queries<E,C>(assembly: &Assembly<E,C>, domain_size: usize) -> Result<(usize, usize)>
where
    E: Engine,
    C: Circuit<E>
{
    let kp_num: _ = &assembly.key_pair.as_ref()?.num;

    let g1_query: usize = {
        let h_query: _ = domain_size;
        let icl_query: _ = kp_num.inputs + kp_num.aux;
        let a_query: _ = kp_num.inputs + kp_num.aux;
        let b_query: _ = kp_num.inputs + kp_num.aux;
        h_query + icl_query + a_query + b_query
    };
    let g2_query: usize = kp_num.inputs + kp_num.aux;

    Ok((g1_query, g2_query))
}
