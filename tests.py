"""Tests for ReliabilityDiagram package"""

import numpy as np
import ReliabilityDiagram as rd

#####################################################################################

def test_dimension_mismatch():
    # test ReliabilityDiagram for dimension mismatch
    nobs = 100   # number of events
    nclim = 30    # number of years considered in climatology
    nfc = 50      # ensemble size of the forecasts

    clima = np.random.rand(nobs,nclim)
    obs = np.random.rand(nobs)
    fcast = np.random.rand(nobs,nfc)

    # First dimension mismatch
    obs1 = np.random.rand(nobs+1)
    fcast1 = np.random.rand(nobs+2,nfc)
    clima1 = np.random.rand(nobs+2,nfc)

    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs1, fcast, clima, 0, 1/4)
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast1, clima, 0, 1/4)
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast, clima1, 0, 1/4)

    # Shape specification
    obs2 = np.random.rand(nobs,2)
    fcast2  = np.random.rand(nobs,nfc,2)
    clima2 = np.random.rand(nobs,nclim,1)

    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs2, fcast, clima, 0, 1/4)
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast2, clima, 0, 1/4)
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast, clima2, 0, 1/4)

    # Weights shape
    weights = np.random.rand(nobs,nfc+1)
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast, clima, 0, 1/4, weights=weights)


def test_nbins():
    # test number of bins
    nobs = 100   # number of events
    nclim = 30    # number of years considered in climatology
    nfc = 50      # ensemble size of the forecasts

    clima = np.random.rand(nobs,nclim)
    obs = np.random.rand(nobs)
    fcast = np.random.rand(nobs,nfc)

    # test nbins type
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast, clima, 0, 1/4,nbins=5.5)
    #test shape bins
    data = rd.ReliabilityDiagram(obs, fcast, clima, 0, 1/4,nbins=10)
    np.testing.assert_array_equal(len(data.bins),10)


def test_bounds_validity():
    # test validity of the bounds
    nobs = 100   # number of events
    nclim = 30    # number of years considered in climatology
    nfc = 50      # ensemble size of the forecasts

    clima = np.random.rand(nobs,nclim)
    obs = np.random.rand(nobs)
    fcast = np.random.rand(nobs,nfc)

    # Bounds outside of [0,1]
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast, clima, -1/4, 1/4)
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast, clima, 3/4, 2)
    # Lower bound > upper bound
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast, clima, 1/4, 0)
    # Special case: lower bound = 0 and upper bound = 1
    data = rd.ReliabilityDiagram(obs, fcast, clima, 0, 1)
    expected_table = np.zeros((5,2))
    expected_table[-1,0] = nobs
    np.testing.assert_array_equal(data.contingency_table(),expected_table)


def test_for_NaN():
    # test NaN values caught
    nobs = 100   # number of events
    nclim = 30    # number of years considered in climatology
    nfc = 50      # ensemble size of the forecasts

    clima = np.random.rand(nobs,nclim)
    obs = np.random.rand(nobs)
    fcast = np.random.rand(nobs,nfc)

    # NaN in observations
    obs1 = obs.copy()
    obs1[0] = np.nan
    clima1 = clima.copy()
    clima1[0] = np.nan
    fcast1 = fcast.copy()
    fcast1[0] = np.nan

    #
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs1, fcast, clima, 0, 1/4)
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast, clima1, 0, 1/4)
    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast1, clima, 0, 1/4)


def test_deterministic_forecast():
    # test deterministic forecasts
    nobs = 100  # number of events
    nclim = 30  # number of years considered in climatology
    nfc = 1    # ensemble size of the forecasts

    clima = np.random.rand(nobs,nclim)
    fcast = np.random.rand(nobs,nfc)
    obs = np.random.rand(nobs)

    np.testing.assert_raises(ValueError, rd.ReliabilityDiagram, obs, fcast, clima, 3/4, 1)


def test_contingency_table():
    # easy check of the contingency table
    nobs = 100  # number of events
    nclim = 30  # number of years considered in climatology
    nfc = 50    # ensemble size of the forecasts
    nbins = 5

    clima = np.random.rand(nobs,nclim)
    fcast = np.random.rand(nobs,nfc)
    obs = np.random.rand(nobs)

    data = rd.ReliabilityDiagram(obs, fcast, clima, 0, 1/4, nbins=nbins)
    # test shape
    np.testing.assert_array_equal(data.contingency_table().shape,np.array([nbins,2]))
    # test sum to number of event
    np.testing.assert_equal(np.sum(data.contingency_table()),nobs)


def test_weights():
    # test weights
    nobs = 100   # number of events
    nclim = 30    # number of years considered in climatology
    nfc = 50      # ensemble size of the forecasts

    clima = np.random.rand(nobs,nclim)
    obs = np.random.rand(nobs)

    # Test 1: equal weights
    fcast = np.random.rand(nobs,nfc)
    weights = np.ones(fcast.shape)
    data = rd.ReliabilityDiagram(obs, fcast, clima, 0, 1/4)
    data_with_weights = rd.ReliabilityDiagram(obs, fcast, clima, 0, 1/4, weights=weights)
    np.testing.assert_array_equal(data.contingency_table(),data_with_weights.contingency_table())

    # Test 2: member with weight 0
    weights = np.ones(fcast.shape)
    weights[:,-1] = 0
    data = rd.ReliabilityDiagram(obs, fcast[:,:-1], clima, 0, 1/4)
    data_with_weights = rd.ReliabilityDiagram(obs, fcast, clima, 0, 1/4, weights=weights)
    np.testing.assert_array_equal(data.contingency_table(),data_with_weights.contingency_table())


def test_closed_ends():
    # test closed_ends option
    nobs = 100   # number of events
    nclim = 30    # number of years considered in climatology
    nfc = 50      # ensemble size of the forecasts

    clima = np.random.rand(nobs,nclim)
    fcast = np.random.rand(nobs,nfc)
    obs = np.random.rand(nobs)
    _, u = np.percentile(clima,[0,100*1/4],axis=1)
    l, _ = np.percentile(np.sort(clima,axis=1),[100*3/4,0],axis=1)
    obs_u = u
    obs_l = l
    fcast_u = np.tile(u,(nfc,1)).T
    fcast_l = np.tile(l,(nfc,1)).T

    # Check upper bound in observations
    data = rd.ReliabilityDiagram(obs_u, fcast, clima, 0, 1/4,closed_ends="both")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=0),np.array([nobs,0]))

    data = rd.ReliabilityDiagram(obs_u, fcast, clima, 0, 1/4,closed_ends="none")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=0),np.array([0,nobs]))

    data = rd.ReliabilityDiagram(obs_u, fcast, clima, 0, 1/4,closed_ends="left")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=0),np.array([0,nobs]))

    data = rd.ReliabilityDiagram(obs_u, fcast, clima, 0, 1/4,closed_ends="right")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=0),np.array([nobs,0]))

    # Check lower bound in observations
    data = rd.ReliabilityDiagram(obs_l, fcast, clima, 3/4, 1,closed_ends="both")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=0),np.array([nobs,0]))

    data = rd.ReliabilityDiagram(obs_l, fcast, clima, 3/4, 1,closed_ends="none")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=0),np.array([0,nobs]))

    data = rd.ReliabilityDiagram(obs_l, fcast, clima, 3/4, 1,closed_ends="left")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=0),np.array([nobs,0]))

    data = rd.ReliabilityDiagram(obs_l, fcast, clima, 3/4, 1,closed_ends="right")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=0),np.array([0,nobs]))

    # Check upper bound in forecasts
    data = rd.ReliabilityDiagram(obs, fcast_u, clima, 0, 1/4,closed_ends="both")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=1),np.array([0,0,0,0,nobs]))

    data = rd.ReliabilityDiagram(obs, fcast_u, clima, 0, 1/4,closed_ends="none")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=1),np.array([nobs,0,0,0,0]))

    data = rd.ReliabilityDiagram(obs, fcast_u, clima, 0, 1/4,closed_ends="left")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=1),np.array([nobs,0,0,0,0]))

    data = rd.ReliabilityDiagram(obs, fcast_u, clima, 0, 1/4,closed_ends="right")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=1),np.array([0,0,0,0,nobs]))

    # Check lower bound in forecasts
    data = rd.ReliabilityDiagram(obs, fcast_l, clima, 3/4, 1,closed_ends="both")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=1),np.array([0,0,0,0,nobs]))

    data = rd.ReliabilityDiagram(obs, fcast_l, clima, 3/4, 1,closed_ends="none")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=1),np.array([nobs,0,0,0,0]))

    data = rd.ReliabilityDiagram(obs, fcast_l, clima, 3/4, 1,closed_ends="left")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=1),np.array([0,0,0,0,nobs]))

    data = rd.ReliabilityDiagram(obs, fcast_l, clima, 3/4, 1,closed_ends="right")
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=1),np.array([nobs,0,0,0,0]))


def test_for_perfect_reliability():
    # test functions for a case with perfect reliability
    nobs = 50   # number of events
    nclim = 30    # number of years considered in climatology
    nfc = 10      # ensemble size of the forecasts

    clima = np.random.rand(nobs,nclim)
    fcast = np.zeros((nobs,nfc))
    fcast[:,0] = 1
    fcast[10:,2:4] = 1
    fcast[20:,4:6] = 1
    fcast[30:,6:8] = 1
    fcast[40:,8:10] = 1
    obs = np.zeros(nobs)
    obs[0] = 1
    obs[10:13] = 1
    obs[20:25] = 1
    obs[30:37] = 1
    obs[40:49] = 1

    data = rd.ReliabilityDiagram(obs, fcast, clima, 3/4, 1,closed_ends="both",nbins=5)

    # Test contingency_table
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=1),10*np.ones(len(data.bins)))
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=0),np.array([np.sum(obs),nobs-np.sum(obs)]))

    # Test observed_frequency
    np.testing.assert_allclose(data.observed_frequency(), data.bins + (data.bins[1]-data.bins[0])/2)

    # Test forecast_attributes (Brier score and reliability)
    np.testing.assert_allclose(data.forecast_attributes()[1],0,atol=1e-16)
    np.testing.assert_allclose(data.forecast_attributes()[0],np.mean((np.sum(fcast,axis=1)/nfc-obs)**2))


def test_for_no_resolution():
    # test functions for a case with no resolution
    nobs = 50   # number of events
    nclim = 30    # number of years considered in climatology
    nfc = 10      # ensemble size of the forecasts


    clima = np.random.rand(nobs,nclim)
    fcast = np.zeros((nobs,nfc))
    fcast[5:,0] = 1
    fcast[10:,1] = 1
    fcast[15:,2] = 1
    fcast[20:,3] = 1
    obs = np.zeros(nobs)
    obs[0::5] = 1

    clima = np.random.rand(nobs,nclim)
    fcast = np.zeros((nobs,nfc))
    fcast[:,0] = 1
    fcast[10:,2:4] = 1
    fcast[20:,4:6] = 1
    fcast[30:,6:8] = 1
    fcast[40:,8:10] = 1
    obs = np.zeros(nobs)
    obs[0::5] = 1

    data = rd.ReliabilityDiagram(obs, fcast, clima, 0.8, 1,closed_ends="both",nbins=5)

    # Test contingency_table
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=1),10*np.ones(len(data.bins)))
    np.testing.assert_array_equal(np.sum(data.contingency_table(),axis=0),np.array([np.sum(obs),nobs-np.sum(obs)]))

    # Test observed_frequency
    np.testing.assert_allclose(data.observed_frequency(),np.repeat(0.2,len(data.bins)))

    # Test forecast_attributes (Brier score and resolution)
    np.testing.assert_allclose(data.forecast_attributes()[2],0,atol=1e-16)
    np.testing.assert_allclose(data.forecast_attributes()[0],np.mean((np.sum(fcast,axis=1)/nfc-obs)**2))


#####################################################################################
# Tests
test_dimension_mismatch()
test_nbins()
test_bounds_validity()
test_for_NaN()
test_deterministic_forecast()
test_contingency_table()
test_weights()
test_closed_ends()
#test_for_perfect_reliability()
test_for_no_resolution()
