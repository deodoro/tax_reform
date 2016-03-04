# -*- coding: utf-8 -*-
"""
This code has been created by Cesare Buiatti, Department of Economics, University of Illinois at Urbana-Champaign, from an original version by Joao B. Duarte.
"""
import tornado.ioloop
import tornado.web
import numpy as np
from scipy.optimize import newton
import random
import json

class Util:
    """
    Utility function class with parameters alpha
    """
    
    def __init__(self, alpha):
        self.alpha = alpha
        
    def total(self, c, h):
        return np.log(c)+self.alpha*np.log(1-h)
    
    def muc(self, c):
        return 1/c

    def muh(self, h):
        return -self.alpha/(1-h)
    
    def mrsctoh(self, c, h):
        return self.muc(h)/self.muh(c)
    
    def mrshtoc(self, c, h):
        return self.muh(h)/self.muc(c)
    
class Prod:
    """
    This class provides the simple production function with parameters A and alpha
    """
    
    def __init__(self,A,theta):
        self.A=A
        self.theta=theta
        
    def total(self,k,h):
        return self.A* k**(self.theta) * h**(1-self.theta)
    
    def mph(self,k,h):
        ktoh=k/h
        return (1-self.theta)*self.A*ktoh**self.theta
    
    def mpk(self,k,h):
        ktoh=k/h
        return self.theta*self.A*ktoh**(self.theta-1)
        
class GrowthModel:
    
    def __init__(self, alpha, beta, A, theta, delta, tauh, tauc, tauk, g=0, G=0, T=10000):
        self. alpha = alpha
        self.beta = beta
        self.A = A
        self.theta = theta
        self.delta = delta
        self.tauh = tauh
        self.tauc = tauc
        self.tauk = tauk
        self.g = g
        self.G = G
        self.T = T
        self.util = Util(alpha= self.alpha)
        self.prod = Prod(A = self.A, theta = self.theta)
    
    def ss(self):
        iss=self.beta**(-1)-1
        rss=iss*(1-self.tauk)**(-1)+self.delta
        ktohss=(self.A*self.theta/rss)**(1/(1-self.theta))
        wss=(1-self.theta)*self.A*ktohss**self.theta
        if self.G==0:        
            hss=(wss*(1-self.tauh)/(1+self.tauc))/(self.alpha*(1-self.g)*self.A*ktohss**self.theta-self.alpha*self.delta*ktohss+wss*(1-self.tauh)/(1+self.tauc))
        elif self.g==0:
            hss=(self.alpha*self.G+wss*(1-self.tauh)/(1+self.tauc))/(self.alpha*self.A*ktohss**self.theta-self.alpha*self.delta*ktohss+wss*(1-self.tauh)/(1+self.tauc))
        kss=ktohss*hss
        css=(1-self.tauh)*wss*(1-hss)/(self.alpha*(1+self.tauc))
        return css,hss,kss,wss,rss,iss
    
    def TaxRev(self):
        css,hss,kss,wss,rss,iss=self.ss()
        TR=self.tauk*(rss-self.delta)*kss+self.tauh*wss*hss+self.tauc*css
        if self.G==0:
            G=self.g*self.prod.total(kss,hss)
        elif self.g==0:
            G=self.G
        T=TR-G
        return np.array([TR,G,T])
    
    def welfare(self):
        W = np.zeros(self.T) # the bigger the number here the more time it takes but the more accurate the sum is 
        for i in range(self.T):
            W[i] = self.beta**i * self.util.total(self.ss()[0], self.ss()[1])
        return sum(W)
        
# Use parameters inserted by the user to estimate the benchamerk model

def bench_economy(params):
    BenchEconomy=GrowthModel(params['alpha'], params['beta'], params['A'], params['theta'], params['delta'], params['tauh'], params['tauc'], params['tauk'], g=params['g']) 

    BenchTaxRev=BenchEconomy.TaxRev()
    return {'consumption': round(BenchEconomy.ss()[0],3),
               'labor': round(BenchEconomy.ss()[1],3),
               'capital': round(BenchEconomy.ss()[2],3),
               'output': round(BenchEconomy.prod.total(BenchEconomy.ss()[2],BenchEconomy.ss()[1]),3),
               'wage': round(BenchEconomy.ss()[3],3),
               'real_return_to_k': round(BenchEconomy.ss()[4],3),
               'real_i_rate': round(BenchEconomy.ss()[5],3),
               'total_revenue': round(BenchTaxRev[0],3),
               'used_for_g_expenditure': round(BenchTaxRev[1],3),
               'used_for_transfers': round(BenchTaxRev[2],3),
               'welfare': round(BenchEconomy.welfare(),3),
               'tauc': round(params['tauc'],3) }

def a_run(params):
    # Find the tax on consumption that holds revenue neutrality
    def objFun(tauc):
#        econ=GrowthModel(alpha,beta,A, theta, delta, tauh, tauc, tauk,G=BenchTaxRev[1])
#        econ=GrowthModel(alpha,beta,A, theta, delta, tauh, tauc, tauk,G=g)
        econ=GrowthModel(params['alpha'], params['beta'], params['A'], params['theta'], params['delta'], params['tauh'], tauc, params['tauk'], G=BenchTaxRev[1])
        return econ.TaxRev()[2]-BenchTaxRev[2]
        
    tauc=newton(objFun,0.5)

    # G does not change in original
    ReformEconomy=GrowthModel(params['alpha'], params['beta'], params['A'], params['theta'], params['delta'], params['tauh'], tauc, params['tauk'], g=params['g']) 
    ReformTaxRev=ReformEconomy.TaxRev()

    return {'consumption': round(ReformEconomy.ss()[0],3),
               'labor': round(ReformEconomy.ss()[1],3),
               'capital': round(ReformEconomy.ss()[2],3),
               'output': round(ReformEconomy.prod.total(ReformEconomy.ss()[2],ReformEconomy.ss()[1]),3),
               'wage': round(ReformEconomy.ss()[3],3),
               'real_return_to_k': round(ReformEconomy.ss()[4],3),
               'real_i_rate': round(ReformEconomy.ss()[5],3),
               'total_revenue': round(ReformTaxRev[0],3),
               'used_for_g_expenditure': round(ReformTaxRev[1],3),
               'used_for_transfers': round(ReformTaxRev[2],3),
               'welfare': round(ReformEconomy.welfare(),3),
               'tauc': round(tauc,3) }


def calculate_params(params):
    for i in params.keys():
        params[i] = float(params[i])
    params['theta']=params['rKtoY']
    params['delta']=params['XtoY']/params['KtoY']
    r=params['rKtoY']/params['KtoY']
    i=(1-params['tauk'])*(r-params['delta'])
    params['beta']=1/(1+i)
    h=float(params['H']/100)
    params['alpha']=((1-h)/h)*(1-params['rKtoY'])*((1-params['tauh'])/(1+params['tauc']))*(1-params['XtoY']-params['GtoY'])**(-1)
    params['A']=1    
    params['g']=params['GtoY']
    return params

#Calibrate
class Calibrate(tornado.web.RequestHandler):
    def post(self):
        print self.request.body
        self.set_header("Content-Type", "application/json")
        params = calculate_params(json.loads(self.request.body))
        self.write(params)

#ARun
class ARun(tornado.web.RequestHandler):
    def post(self):
        print self.request.body
        self.set_header("Content-Type", "application/json")
        params = calculate_params(json.loads(self.request.body))
        self.write(a_run(params))

# TEST
class Benchmark(tornado.web.RequestHandler):
    def post(self):
        print self.request.body
        self.set_header("Content-Type", "application/json")
        params = calculate_params(json.loads(self.request.body))
        self.write(bench_economy(params))

class MainHandler(tornado.web.RequestHandler):
    def get(self):
        print 'test'
        self.render('static/index.html')

if __name__ == "__main__":
    settings = {
        "debug": True,
        "static_path": "static"
    }

    application = tornado.web.Application([(r'/bench', Benchmark),
                                           (r'/run', ARun),
                                           (r'/calibrate', Calibrate),
                                           (r'/', MainHandler)],
                                           **settings)
    application.listen(8080)
    tornado.ioloop.IOLoop.instance().start()
