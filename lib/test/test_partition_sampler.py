from math import log, exp
from random import betavariate as beta_rvs

from pyclone.sampler import Partition, sample_alpha, sample_partition
from pyclone.utils import binomial_rvs, gamma_rvs, log_binomial_pdf, log_space_normalise, log_beta_pdf, discrete_rvs, \
    log_binomial_coefficient, log_beta

from pyclone.test.crp import sample_from_crp
from pyclone.test.geweke import compare_trace, mean

class BetaBaseMeasure(object):
    def __init__(self, a, b):
        self.a = a
        self.b = b
        
    def random(self):
        return beta_rvs(self.a, self.b)

class BinomialDensity(object):    
    def log_p(self, data, p):
        x = data[0]
        n = data[1]
        
        return log_binomial_pdf(x, n, p)

def log_beta_binomial_pdf(k, n, a, b):
    return log_binomial_coefficient(n, k) + log_beta(a + k, b + n - k) - log_beta(a, b)

#=======================================================================================================================
# 
#=======================================================================================================================
def sample_cell_values(data, partition, base_measure):
    for cell in partition.cells:
        a = base_measure.a
        b = base_measure.b
        
        for item in cell.items:
            a += data[item][0]
            b += data[item][1] - data[item][0]
        
        cell.value = beta_rvs(a, b)

def sample_partition_conjugate(alpha, data, old_partition, base_mesure, cluster_density):
    partition = old_partition.copy()
    
    for item, data_point in enumerate(data):
        old_cell_index = partition.labels[item]
        
        partition.remove_item(item, old_cell_index)
        
        partition.remove_empty_cells()
        
        log_p = []    
        
        for cell in partition.cells:
            cluster_log_p = cluster_density.log_p(data_point, cell.value)
            
            counts = cell.size
            
            log_p.append(log(counts) + cluster_log_p)

        cluster_log_p = log_beta_binomial_pdf(data_point[0], data_point[1], base_mesure.a, base_mesure.b)
        
        log_p.append(log(alpha) + cluster_log_p)
        
        log_p = log_space_normalise(log_p)    
        
        p = [exp(x) for x in log_p]
        
        new_cell_index = discrete_rvs(p)
        
        if new_cell_index == partition.number_of_cells:
            partition.add_cell(base_mesure.random())
            
        partition.add_item(item, new_cell_index)
    
    return partition

#=======================================================================================================================
# 
#=======================================================================================================================
def draw_parameters(gamma_a, gamma_b, base_measure, size): 
    partition = Partition()
    
    alpha = gamma_rvs(gamma_a, gamma_b)
    
    labels, values = sample_from_crp(alpha, size, base_measure)
    
    for k, value in enumerate(values):
        partition.add_cell(value)
        
        for item, cell_index in enumerate(labels):
            if cell_index == k:
                partition.add_item(item, cell_index)

    return alpha, partition

def draw_data(partition, n):
    data = []
    
    for p in partition.item_values:
        data.append((binomial_rvs(n, p), n))
    
    return data

#=======================================================================================================================
# 
#=======================================================================================================================
def marginal_conditional_simulator(gamma_a, gamma_b, beta_a, beta_b, size, num_iters): 
    params = []
    
    base_measure = BetaBaseMeasure(beta_a, beta_b)
    
    for _ in range(num_iters):
        alpha, partition = draw_parameters(gamma_a, gamma_b, base_measure, size)
        
        params.append({'alpha' : alpha, 'p' : partition.item_values})
    
    return params

def successive_conditional_simulator(gamma_a, gamma_b, beta_a, beta_b, n, size, num_iters):
    params_sample = []
    
    base_measure = BetaBaseMeasure(beta_a, beta_b)
    
    cluster_density = BinomialDensity()
    
    alpha, partition = draw_parameters(gamma_a, gamma_b, base_measure, size)
    
    for _ in range(num_iters):
        data = draw_data(partition, n)
        
        new_alpha = sample_alpha(alpha,
                                 partition.number_of_cells,
                                 partition.number_of_items,
                                 a=gamma_a,
                                 b=gamma_b)        
        
        new_partition = sample_partition_conjugate(new_alpha, data, partition, base_measure, cluster_density)
               
        sample_cell_values(data, new_partition, base_measure)
        
        params = {'alpha' : new_alpha, 'p' : new_partition.item_values}
        
        params_sample.append(params)
        
        alpha = new_alpha
        
        partition = new_partition
    
    return params_sample

def dp_sampler(data, gamma_a, gamma_b, beta_a, beta_b, num_iters):
    params_sample = []
    
    base_measure = BetaBaseMeasure(beta_a, beta_b)
    
    cluster_density = BinomialDensity()
    
    partition = Partition()
    
    alpha = gamma_rvs(gamma_a, gamma_b)
    
    for item, _ in enumerate(data):
        partition.add_cell(base_measure.random())
        partition.add_item(item, item)
    
    for _ in range(num_iters):
        alpha = sample_alpha(alpha,
                             partition.number_of_cells,
                             partition.number_of_items,
                             a=gamma_a,
                             b=gamma_b)        
        
        partition = sample_partition_conjugate(alpha, data, partition, base_measure, cluster_density)
               
        sample_cell_values(data, partition, base_measure)
        
        params = {'alpha' : alpha, 'p' : partition.item_values, 'labels' : partition.labels}
        
        params_sample.append(params)
    
    return params_sample

#=======================================================================================================================
# 
#=======================================================================================================================
if __name__ == "__main__":
    n = 100
    
    gamma_a = 1
    gamma_b = 1
    
    beta_a = 1
    beta_b = 1
    
    size = 10
    num_iters = 1000000
    
    burnin = 1
    thin = 100
    
#    params_1 = marginal_conditional_simulator(gamma_a, gamma_b, beta_a, beta_b, size, num_iters)
#    params_2 = marginal_conditional_simulator(gamma_a, gamma_b, beta_a, beta_b, size, num_iters)
    
    params_1 = marginal_conditional_simulator(gamma_a, gamma_b, beta_a, beta_b, size, num_iters)
    params_2 = successive_conditional_simulator(gamma_a, gamma_b, beta_a, beta_b, n, size, num_iters)    
    
    trace_1 = [x['alpha'] for x in params_1][burnin::thin]
    trace_2 = [x['alpha'] for x in params_2][burnin::thin]
    
    print compare_trace(trace_1, trace_2, lambda x: x)
    print compare_trace(trace_1, trace_2, lambda x: x ** 2)
    
    for i in range(size):
        trace_1 = [x['p'][i] for x in params_1][burnin::thin]
        trace_2 = [x['p'][i] for x in params_2][burnin::thin]
        
        print compare_trace(trace_1, trace_2, lambda x: x)
        print compare_trace(trace_1, trace_2, lambda x: x ** 2)

#    data = [
#            (20, 100),
#            (20, 100),
#            (20, 100),
#            (50, 100),
#            (50, 100),
#            (50, 100),
#            (50, 100),
#            (99, 100),
#            (99, 100),
#            ]
#    
#    params = dp_sampler(data, gamma_a, gamma_b, beta_a, beta_b, num_iters)
#    
#    for i in range(len(data)):
#        trace = [x['p'][i] for x in params][burnin::thin]
#        
#        print mean(trace)
#    
#    trace = [len(set(x['labels'])) for x in params[burnin::thin]]
#    
#    print mean(trace)