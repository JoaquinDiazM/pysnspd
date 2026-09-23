import numpy as np
import pytest
from pysnspd.experimental.energy_catalog import energy_at_count_batch
from pysnspd.experimental.seeded_count_catalog import invert_count_seeded


@pytest.mark.parametrize('seed_kind',['near','reversed','outside'])
def test_same_causal_inverse_from_arbitrary_trial_guesses(seed_kind):
    x=np.r_[np.geomspace(1e-8,.01,7),np.linspace(.02,3,9)]
    exact=energy_at_count_batch(x,delta=.995,gamma=.0055,eta=1e-8)
    seed={'near':exact*1.00001,'reversed':exact[::-1],'outside':np.full(len(x),-1.)}[seed_kind]
    corrected,diag=invert_count_seeded(x,.995,.0055,1e-8,seed)
    np.testing.assert_allclose(corrected,exact,rtol=1e-10,atol=2e-12)
    assert np.all(np.diff(corrected)>0) and diag['spectrum_points']>0
    if seed_kind=='outside':assert diag['seed_bracket_replacements']==len(x)


def test_invalid_inputs_rejected():
    with pytest.raises(ValueError):invert_count_seeded([0.,1.],1.,.01,1e-8,[1.,2.])
    with pytest.raises(ValueError):invert_count_seeded([.1,1.],1.,0.,1e-8,[1.,2.])
