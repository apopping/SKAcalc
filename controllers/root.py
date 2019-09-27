# -*- coding: utf-8 -*-
"""Main Controller"""
from tg import expose, flash, require, url, lurl
from tg import request, redirect, tmpl_context
from tg.i18n import ugettext as _, lazy_ugettext as l_
from tg.exceptions import HTTPFound
from tg import predicates
from ska_calculator import model
from ska_calculator.controllers.secure import SecureController
from ska_calculator.model import DBSession
from tgext.admin.tgadminconfig import BootstrapTGAdminConfig as TGAdminConfig
from tgext.admin.controller import AdminController

from ska_calculator.lib.base import BaseController
from ska_calculator.controllers.error import ErrorController

import os
import sys
import time

import calculate_sensitivity as cs
from read_parameters import read_parameters

__all__ = ['RootController']




class RootController(BaseController):
    """
    The root controller for the ska_calculator application.

    All the other controllers and WSGI applications should be mounted on this
    controller. For example::

        panel = ControlPanelController()
        another_app = AnotherWSGIApplication()

    Keep in mind that WSGI applications shouldn't be mounted directly: They
    must be wrapped around with :class:`tg.controllers.WSGIAppController`.

    """
    secc = SecureController()
    admin = AdminController(model, DBSession, config_type=TGAdminConfig)

    error = ErrorController()


    #movie = MovieController()
    @expose()
    def index(self):
        return "<h1>Hello World</h1>"

    @expose()
    def _default(self, *args, **kw):
        return "This page is not ready"

    @expose("ska_calculator.templates.input")
    def input2(self):
        return{}


    @expose("ska_calculator.templates.sample")
    def output(self, **kw):
        current_time = time.strftime('%Y%m%d%H%M%S')
        # read the pre-defined ska parameters
        parameter_file = 'ska_calculator/data/ska_parameters.par'
        ska_par = read_parameters(parameter_file)
        # check whether input is valid
        check, message = cs.check_input(kw, ska_par)
        if check == 'Fail':
            return message
        else:
            # generate the uv file
            cs.make_uv(kw,current_time, ska_par)
            # generate the noise image and psf
            cs.make_image(kw,current_time, ska_par)
            # make a png plot of the psf
            # make a directory for the output:
            #public_dir = 'ska_calculator/public/data/' + current_time + '/'
            #os.system('mkdir ' +  public_dir)
            cs.plot_psf(kw,current_time)
            cs.plot_noise(kw,current_time)
            # plot the visibilities:
            # the fucntion below is disabled, it is working, however miriad has a memory limit which becomes problematic for LOW where the number of baselines exceeds teh limit to plot every point
            
            #cs.plot_uv(kw,current_time, ska_par)
            # generate some results
            #current_time4 = time.strftime('%Y%m%d%H%M%S')
            mydata = cs.ska_sens(kw, current_time)
            # clean the temporary data:
            os.system('rm -rf /tmp/ska_calculator/data/runs/' + current_time)
            return mydata
