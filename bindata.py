# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#import sqlite3 as sql
import numpy as np
import astropy as ap
import time
import sys
import scipy
import scipy.stats
import matplotlib.pyplot as plt
import mysql.connector as sql
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import pandas
import pandas.io.sql as pandasql
import seaborn
import matplotlib
from sqlalchemy import and_, or_

#s_root = "ngc5548_1e00_obs_2"
s_user = "root"
s_password = "password"

s_file = "test_orig"

### TRY OPENING THE DATABASE ###
print ("Opening database '{}'...".format(s_file))

db_engine = None
try:
    db_engine = sqlalchemy.create_engine("sqlite:///{}.db".format(s_file))
except sqlalchemy.exc.SQLAlchemyError as e:
    print(e)
    sys.exit(1)

### DOES IT ALREADY EXIST? ###
print ("Searching for table 'Photons'...")
Base = sqlalchemy.ext.declarative.declarative_base()
class Spectrum(Base):
    __tablename__ = "Spectra"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    angle = sqlalchemy.Column(sqlalchemy.Float)

class Origin(Base):
    __tablename__ = "Origins"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String)

class Photon(Base):
    __tablename__ = "Photons"
    id = sqlalchemy.Column(sqlalchemy.Integer, primary_key=True, autoincrement=True)
    Wavelength = sqlalchemy.Column(sqlalchemy.Float)
    Weight = sqlalchemy.Column(sqlalchemy.Float)
    X = sqlalchemy.Column(sqlalchemy.Float)
    Y = sqlalchemy.Column(sqlalchemy.Float)
    Z = sqlalchemy.Column(sqlalchemy.Float)
    ContinuumScatters = sqlalchemy.Column(sqlalchemy.Integer)
    ResonantScatters = sqlalchemy.Column(sqlalchemy.Integer)
    Delay = sqlalchemy.Column(sqlalchemy.Float)
    Spectrum = sqlalchemy.Column(sqlalchemy.Integer)
    Origin = sqlalchemy.Column(sqlalchemy.Integer)
    Resonance = sqlalchemy.Column(sqlalchemy.Integer)
    Origin_matom = sqlalchemy.Column(sqlalchemy.Boolean)

Session = sqlalchemy.orm.sessionmaker(bind=db_engine)
dbc = Session()

start = time.clock()

try:
    dbc.query(Photon).first()
    # If so, we go with what we've found.
    print("Found existing filled photon database '{}'".format(s_file))
except sqlalchemy.exc.SQLAlchemyError as e:
    print(e)
    # If not, we populate from the delay dump file. This bit is legacy!
    print("No existing filled photon database, reading from file '{}.delay_dump'".format(s_file))
    Base.metadata.create_all(db_engine)

    delay_dump = open("{}.delay_dump".format(s_file), 'r')
    for line in delay_dump:
        if line.startswith('#'): 
            continue
        values = line.split()
        if len(values) < 12:
            continue
        
        for index,value in enumerate(values):
            values[index] = float(value)
        del values[0]
        del values[8]
        values[5] = (values[5] - values[6])
        matom_bool = False
        if(values[9]>= 10):
            values[9] = values[9] - 10
            matom_bool = True

        dbc.add(Photon(Wavelength=values[0], Weight=values[1], X=values[2], Y=values[3], Z=values[4],
                        ContinuumScatters=values[5], ResonantScatters=values[6], Delay=values[7],
                        Spectrum=values[8], Origin=values[9], Resonance=values[10], Origin_matom = matom_bool))
    dbc.commit()

print("Input took {:.2f} seconds".format(time.clock()-start))

query = dbc.query(Photon.Wavelength, Photon.Delay, Photon.Weight)

# HOW TO BUILD A QUERY:
query = query.filter(and_(Photon.Spectrum == 0))
#query = query.filter(and_(Photon.Wavelength > 1475, Photon.Wavelength < 1625))
#query = query.filter(and_(or_(Photon.ResonantScatters > 0, Photon.Origin == 3, Photon.Origin_matom == True)))
#query = query.filter(and_(Photon.Resonance.in_([380, 381, 382, 416])))


start = time.clock()
data = pandasql.read_sql_query(query.subquery(), db_engine)
seconds_to_days = 60 * 60 * 24
data['Delay'] = data['Delay'].apply(lambda x: x/seconds_to_days)
print("PANDAS: {}".format(time.clock() - start))
seaborn.set_style("white")
seaborn.set_style("ticks")
seaborn.set_context("paper")

def weighted_hex(x,y, **kwargs):
    pandas.DataFrame.plot.hexbin(x,y,C=data.Weight, **kwargs)

jaxis = seaborn.JointGrid("Wavelength", "Weight", data=data)
jaxis.plot_joint(weighted_hex)
#jaxis = seaborn.jointplot("Wavelength", "Delay", data=data, space=0, kind="hex", stat_func=None, joint_kws=dict(C="Weight"), shade_lowest=False)
#caxis = jaxis.fig.add_axes([1,.34, .01, .2])
#seaborn.plt.ylim(ymin=0)
#seaborn.plt.colorbar(cax=caxis)
jaxis.set_axis_labels(r"Wavelength ($\AA$)", r"Delay $\tau$ (days)")
seaborn.plt.savefig("{}.eps".format(s_file))
