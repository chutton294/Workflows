/*
Copyright (C) 2015 Christopher Hutton.
Copyright (C) 2010-2013 Jon Herman, Josh Kollat, and others.

This file is part of HymodR.

HymodR is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
HymodR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License
along with Hymod.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <Rcpp.h>

#include <hymodR.h>

using namespace Rcpp;

struct HyMod hymod;

int nDays;

// [[Rcpp::export]]
NumericVector Rhymod(NumericVector rain, NumericVector evap, NumericVector temp, NumericVector params){
  
    nDays = rain.size();
    //assign input data to the appropriate vectors
    hymod.data.nDays = nDays;
    hymod.parameters.Nq = 3; // number of quickflow reservoirs
    hymod.parameters.Kv = 1.0; // vegetation parameter
    hymod_allocate(nDays); //allocate size of vectors to store input data
    
    //assign incoming numeric vectors to the appropriate hymod.data pointers
    //hymod.data.precip = new double[nDays];
    for (int i = 0; i < nDays; i++){
        hymod.data.precip[i] = rain[i];
        hymod.data.avgTemp[i] = temp[i];
        hymod.fluxes.PE[i] = evap[i];
    }
    
    //assign incoming parameters from R to model parameters
    int nPar = params.size();
    double* parameters = new double [nPar]; 
    for (int i = 0; i < nPar; i++){
        parameters[i] = params[i];
    }
    
    // Run model with this set of parameters
    calc_hymod(parameters);
       
    //convert output discharge to a numeric vector for output back into R
    int tt = 4 + hymod.parameters.Nq;
    tt = 1;
    NumericMatrix store(nDays, tt); 
    for (int i = 0; i < nDays; i++){
        store(i,0) = hymod.fluxes.Q[i];
        //store(i,1) = hymod.fluxes.Qq[i];
        //store(i,3) = hymod.fluxes.OV[i];
        //store(i,2) = hymod.fluxes.Qs[i];
        //store(i,5) = hymod.fluxes.snow[i];
        //store(i,6) = hymod.fluxes.melt[i];
        //store(i,7) = hymod.fluxes.PE[i];
        //store(i,8) = hymod.states.XHuz[i];
        //store(i,9) = hymod.states.XCuz[i];
        //store(i,10) = hymod.states.Xs[i];
        //store(i,11) = hymod.states.snow_store[i];
        //for(int j = 0; j <  hymod.parameters.Nq;j++){
              //store(i,4+j) = hymod.states.Xq[i][j]; 
        //} 
    } 
             
    hymod_delete(nDays);
    delete[] parameters;
          
   
  return(store);
}

//This is the function that gets called to evaluate each parameter set
void calc_hymod(double* parameters)
{
    // assign parameter values for this run
    // Rate constants Ks and Kq should be specified in units of time-1, 0 < Ks < Kq < 1
    hymod.parameters.Ks    = 1/parameters[0];
    hymod.parameters.Kq    = 1/parameters[1];
    hymod.parameters.DDF   = parameters[2];
    hymod.parameters.Tb    = parameters[3];
    hymod.parameters.Tth   = parameters[4];
    hymod.parameters.alpha = parameters[5];
    hymod.parameters.B     = parameters[6];
    hymod.parameters.Huz   = parameters[7];
    hymod.parameters.Cpar = hymod.parameters.Huz / (1.0 + hymod.parameters.B); // max capacity of soil moisture tank

    // Reinitialize everything to zero
    zero_states_and_fluxes(nDays);

    //Run Model for Simulation Period
    int dataDay;
    for (int modelDay = 0; modelDay < nDays; modelDay++)
    {
        //Since used as an index, we need to convert to zero indexing
        //dataDay = startingIndex + modelDay;
        dataDay = modelDay;

        // Run snow model to find effective precip for this timestep
        hymod.fluxes.effPrecip[modelDay] = snowDD(modelDay, dataDay);

        // Run Pdm soil moisture accounting including evapotranspiration
        PDM_soil_moisture(modelDay, dataDay);

        // Run Nash Cascade routing of quickflow component
        double new_quickflow = hymod.parameters.alpha * hymod.fluxes.OV[modelDay];
        hymod.fluxes.Qq[modelDay] = Nash(hymod.parameters.Kq, hymod.parameters.Nq, new_quickflow, hymod.states.Xq[modelDay]);

        // Run Nash Cascade routing of slowflow component
        double new_slowflow = (1.0-hymod.parameters.alpha) * hymod.fluxes.OV[modelDay];
        hymod.fluxes.Qs[modelDay] = Nash(hymod.parameters.Ks, 1, new_slowflow, &hymod.states.Xs[modelDay]);

        // Set the intial states of the next time step to those of the current time step
        if (modelDay < nDays-1)
        {
            hymod.states.XHuz[modelDay+1] = hymod.states.XHuz[modelDay];
            hymod.states.Xs[modelDay+1]   = hymod.states.Xs[modelDay];
            hymod.states.snow_store[modelDay+1] = hymod.states.snow_store[modelDay];

            for(int m = 0; m < hymod.parameters.Nq; m++)
                hymod.states.Xq[modelDay+1][m] = hymod.states.Xq[modelDay][m];
        }

        hymod.fluxes.Q[modelDay] = hymod.fluxes.Qq[modelDay] + hymod.fluxes.Qs[modelDay];
    }

    return;
}


void zero_states_and_fluxes(int ndays)
{
    for (int k=0; k < ndays; k++)
    {
        hymod.fluxes.snow[k]  = 0.0;
        hymod.fluxes.melt[k]  = 0.0;
        hymod.fluxes.effPrecip[k] = 0.0;
        hymod.fluxes.AE[k] = 0.0;
        hymod.fluxes.OV[k] = 0.0;
        hymod.fluxes.Qq[k] = 0.0;
        hymod.fluxes.Qs[k] = 0.0;
        hymod.fluxes.Q[k]  = 0.0;

        hymod.states.snow_store[k] = 0.0;
        hymod.states.XHuz[k]      = 0.0;
        hymod.states.XCuz[k] = 0.0;
        hymod.states.Xs[k] = 0.0;
        for (int m=0; m < hymod.parameters.Nq; m++) hymod.states.Xq[k][m] = 0.0;
    }
    return;
}

// Allocate time series arrays for states, fluxes, and forcing data
void hymod_allocate(int ndays)
{
    hymod.states.snow_store     = new double [ndays];
    hymod.states.XHuz = new double [ndays];
    hymod.states.XCuz = new double [ndays];
    hymod.states.Xs   = new double [ndays];
    hymod.states.Xq   = new double* [ndays];
    for (int i=0; i < ndays; i++) 
        hymod.states.Xq[i] = new double[hymod.parameters.Nq];

    hymod.fluxes.snow           = new double [ndays];
    hymod.fluxes.melt           = new double [ndays];
    hymod.fluxes.effPrecip      = new double [ndays];
    hymod.fluxes.AE             = new double [ndays];
    hymod.fluxes.OV             = new double [ndays];
    hymod.fluxes.Qq             = new double [ndays];
    hymod.fluxes.Qs             = new double [ndays];
    hymod.fluxes.Q              = new double [ndays];
    hymod.fluxes.PE = new double [ndays];
    
    //added to store the input data for the model run.
    hymod.data.precip   = new double[nDays];
    hymod.data.evap     = new double[nDays];
    hymod.data.flow     = new double[nDays];
    hymod.data.maxTemp  = new double[nDays];
    hymod.data.minTemp  = new double[nDays];
    hymod.data.avgTemp  = new double[nDays];    
    
}

// clean up memory
void hymod_delete(int ndays) 
{
    delete[] hymod.fluxes.snow;
    delete[] hymod.fluxes.melt;
    delete[] hymod.states.snow_store;
    
    delete[] hymod.fluxes.effPrecip;
    delete[] hymod.states.XHuz;
    delete[] hymod.states.XCuz;
    delete[] hymod.states.Xs;
    delete[] hymod.fluxes.AE;
    delete[] hymod.fluxes.OV;
    delete[] hymod.fluxes.Qq;
    delete[] hymod.fluxes.Qs;
    delete[] hymod.fluxes.Q;
    
    for (int i=0; i < ndays; i++) delete[] hymod.states.Xq[i];
    delete[] hymod.states.Xq;
    
    delete[] hymod.data.precip;
    delete[] hymod.data.avgTemp;
    delete[] hymod.fluxes.PE;
    delete[] hymod.data.evap;     
    delete[] hymod.data.flow;     
    delete[] hymod.data.maxTemp; 
    delete[] hymod.data.minTemp;      
    
}

void PDM_soil_moisture(int modelDay, int dataDay)
{
    double Cbeg, OV2, PPinf, Hint, Cint, OV1; // temporary variables for intermediate calculations
    
    // Storage contents at begining
    Cbeg = hymod.parameters.Cpar * (1.0 - pow(1.0-(hymod.states.XHuz[modelDay]/hymod.parameters.Huz),1.0+hymod.parameters.B));

    // Compute overflow from soil moisture storage element
    OV2 = max(0.0, hymod.fluxes.effPrecip[modelDay] + hymod.states.XHuz[modelDay] - hymod.parameters.Huz);

    // Remaining net rainfall
    PPinf = hymod.fluxes.effPrecip[modelDay] - OV2;

    // New actual height in the soil moisture storage element
    Hint = min(hymod.parameters.Huz, hymod.states.XHuz[modelDay] + PPinf);

    // New storage content
    Cint = hymod.parameters.Cpar*(1.0-pow(1.0-(Hint/hymod.parameters.Huz),1.0+hymod.parameters.B));

    // Additional effective rainfall produced by overflow from stores smaller than Cmax
    OV1 = max(0.0, PPinf + Cbeg - Cint);

    // Compute total overflow from soil moisture storage element
    hymod.fluxes.OV[modelDay] = OV1 + OV2;
    
    // Compute actual evapotranspiration
    hymod.fluxes.AE[modelDay] = min(Cint, (Cint/hymod.parameters.Cpar)*hymod.fluxes.PE[modelDay]*hymod.parameters.Kv);
    
    // Storage contents and height after ET occurs
    hymod.states.XCuz[modelDay] = max(0.0, Cint - hymod.fluxes.AE[modelDay]);
    hymod.states.XHuz[modelDay] = hymod.parameters.Huz*(1.0-pow(1.0-(hymod.states.XCuz[modelDay]/hymod.parameters.Cpar),1.0/(1.0+hymod.parameters.B)));

    return;

}

double Nash(double K, int N, double Qin, double *X)
{
    //Initialization
    double *OO = new double[N];
    double Qout;                       //Flow out of series of reservoirs
    
    //Loop through reservoirs - outflow cannot be greater than store as long as K < 0.
    for (int Res = 0; Res < N; Res++)
    {
        OO[Res] = K*X[Res];
        X[Res]  = X[Res] - OO[Res];

        if (Res==0) X[Res] = X[Res] + Qin; 
        else        X[Res] = X[Res] + OO[Res-1];
    }

    // The outflow from the cascade is the outflow from the last reservoir
    Qout = OO[N-1];
    delete[] OO;
    return Qout;
}

double snowDD(int modelDay, int dataDay)
{
    double Qout; // effective precip after freezing/melting

    //If temperature is lower than threshold, precip is all snow
    if (hymod.data.avgTemp[dataDay] < hymod.parameters.Tth)
    {
        hymod.fluxes.snow[modelDay] = hymod.data.precip[dataDay];
        Qout = 0.0;
    }
    else //Otherwise, there is no snow and it's all rain
    {
        hymod.fluxes.snow[modelDay] = 0.0;
        Qout = hymod.data.precip[dataDay];
    }

    //Add to the snow storage for this day
    hymod.states.snow_store[modelDay] += hymod.fluxes.snow[modelDay];

    //Snow melt occurs if we are above the base temperature (either a fraction of the store, or the whole thing)
    if (hymod.data.avgTemp[dataDay] > hymod.parameters.Tb)
    {
        hymod.fluxes.melt[modelDay] = min(hymod.parameters.DDF*(hymod.data.avgTemp[dataDay]-hymod.parameters.Tb), hymod.states.snow_store[modelDay]);
    }
    //Otherwise, snowmelt is zero
    else
    {
        hymod.fluxes.melt[modelDay] = 0.0;
    }

    //Update the snow storage depending on melt
    hymod.states.snow_store[modelDay] -= hymod.fluxes.melt[modelDay];
    if(hymod.states.snow_store[modelDay] < 0.0) hymod.states.snow_store[modelDay] = 0.0;

    //Qout is any rain + snow melt
    Qout += hymod.fluxes.melt[modelDay];

    return Qout;
}










