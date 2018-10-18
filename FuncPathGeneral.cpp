/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include <cmath>

#include "Function.h"
#include "ActionRegister.h"

#include <string>
#include <cstring>
#include <iostream>

#include "tools/IFile.h"

using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION FUNCPATHGENERAL
/*
This function calculates path collective variables for any combination of variables. 

Detailed documentation to follow at some point.

*/
//+ENDPLUMEDOC
   
class FuncPathGeneral : public Function {
  double lambda;
  int neigh_size;
  double neigh_stride;
    
  vector<double> coefficients;
  vector< vector<double> > path_cv_values;
  
  // for faster calculation
  vector<double> expdists;
  
  // for calculating derivatives
  vector< vector<double> > numerators;
  vector<double> s_path_ders;
  vector<double> z_path_ders;
    
  // for handling periodicity
  vector<double> domains;
    
  string reference;
  vector<int> columns;
    
  vector< pair<int,double> > neighpair;
  vector <Value*> allArguments; 
    
  // methods
  void loadReference();

struct pairordering {
  bool operator ()(pair<int, double> const& a, pair<int, double> const& b) {
    return (a).second < (b).second;
  }
};

public:
  explicit FuncPathGeneral(const ActionOptions&);
// active methods:
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(FuncPathGeneral,"FUNCPATHGENERAL")

void FuncPathGeneral::loadReference() {
  IFile input;
  input.open(reference);
  if (!input) plumed_merror("Could not open reference file!");
  while (input)
  {
    vector<string> strings;
    Tools::getParsedLine(input, strings);
    if (strings.empty()) continue;
    vector<double> colvarLine;
    double value;
    int max = columns.empty() ? strings.size() : columns.size();
    for (int i = 0; i < max; ++i)
    {
      int col = columns.empty() ? i : columns[i];
      // if no columns entered, ignore the first (time) and take the rest
      if (columns.empty() && i == 0) continue;
      
      Tools::convert(strings[col], value);
      colvarLine.push_back(value);
    }
    path_cv_values.push_back(colvarLine);
  }
}
    
void FuncPathGeneral::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","LAMBDA","the lambda parameter is needed for smoothing, is in the units of plumed");
  keys.add("compulsory","COEFFICIENTS","the coefficients to be assigned to each CV");
  keys.add("compulsory","REFERENCE","the colvar file is needed to provide the various milestones");
  keys.add("optional","COLUMNS","columns in the reference colvar file specifying the CVs");
  keys.add("optional","NEIGH_SIZE","size of the neighbor list");
  keys.add("optional","NEIGH_STRIDE","how often the neighbor list needs to be calculated in time units");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("s","default","the position on the path");
  keys.addOutputComponent("z","default","the distance from the path");
}
    
FuncPathGeneral::FuncPathGeneral(const ActionOptions&ao):
Action(ao),
Function(ao),
neigh_size(-1),
neigh_stride(-1.)
{
  parse("LAMBDA",lambda);
  parse("NEIGH_SIZE",neigh_size);
  parse("NEIGH_STRIDE",neigh_stride);
  parse("REFERENCE",reference);
  parseVector("COEFFICIENTS",coefficients);
  parseVector("COLUMNS",columns);
  checkRead();
  log.printf("  lambda is %f\n",lambda);
  if (getNumberOfArguments() != coefficients.size()) plumed_merror("The numbers of coefficients and CVs are different!");
  if (!columns.empty()) {
    if (columns.size() != coefficients.size()) plumed_merror("The numbers of coefficients and columns are different!");
  }
  log.printf("  Consistency check completed! Your path cvs look good!\n"); 
  
  // load the reference colvar file
  loadReference();
  
  // do some neighbor printout
  if(neigh_stride>0. || neigh_size>0){
           if(neigh_size>path_cv_values.size()){
           log.printf(" List size required ( %d ) is too large: resizing to the maximum number of arg required: %d  \n", neigh_size, getNumberOfArguments());
           neigh_size=path_cv_values.size();
           }
           log.printf("  Neighbor list enabled: \n");
           log.printf("                size   :  %d elements\n",neigh_size);
           log.printf("                stride :  %f time \n",neigh_stride);
  }else{
           log.printf("  Neighbor list NOT enabled \n");
  }

  addComponentWithDerivatives("s"); componentIsNotPeriodic("s");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
    
  // initialise vectors
  vector<double> temp (coefficients.size());
  for (unsigned i = 0; i < path_cv_values.size(); ++i) {
    numerators.push_back(temp);
    expdists.push_back(0.);
    s_path_ders.push_back(0.);
    z_path_ders.push_back(0.);
  }

  // now backup the arguments 
  for(unsigned i=0;i<getNumberOfArguments();i++)allArguments.push_back(getPntrToArgument(i)); 
    
  // get periodic domains, negative for not periodic, stores half the domain length (maximum difference)
  for (unsigned i = 0; i < allArguments.size(); ++i) {
    if (allArguments[i]->isPeriodic()) {
      double min_lim, max_lim;
      allArguments[i]->getDomain(min_lim, max_lim);
      domains.push_back((max_lim - min_lim) / 2);
    }
    else
      domains.push_back(-1.);
  }
}
    
// calculator
void FuncPathGeneral::calculate() {
 // log.printf("NOW CALCULATE! \n");
  double s_path=0.;
  double partition=0.;
  double tmp, value, diff, expdist, s_der, z_der;
  int ii;
    
  typedef vector< pair< int,double> >::iterator pairiter;
    
  for (pairiter it = neighpair.begin(); it != neighpair.end(); ++it) {
    (*it).second = 0.;
  }
    
  if (neighpair.empty()) { // at first step, resize it
    neighpair.resize(path_cv_values.size());
    for(unsigned i=0; i<path_cv_values.size(); ++i) neighpair[i].first=i; 
  }
  
  Value* val_s_path=getPntrToComponent("s");
  Value* val_z_path=getPntrToComponent("z");
    
  for(int j = 0; j < allArguments.size(); ++j) {
    value = allArguments[j]->get();
    for (pairiter it = neighpair.begin(); it != neighpair.end(); ++it) {
      diff = (value - path_cv_values[(*it).first][j]);
      if (domains[j] > 0) {
        if (diff > domains[j]) diff -= 2 * domains[j];
        if (diff < -domains[j]) diff += 2 * domains[j];
      }
      (*it).second += Tools::fastpow(coefficients[j] * diff, 2);
      numerators[(*it).first][j] = 2 * Tools::fastpow(coefficients[j], 2) * diff;
    }
  }
    
  for (pairiter it = neighpair.begin(); it != neighpair.end(); ++it) {
    //(*it).second = sqrt((*it).second);  // changed
    expdist = exp(-lambda*(*it).second);
    expdists[(*it).first] = expdist;
    s_path += ((*it).first + 1) * expdist;
    partition += expdist;
  }
      
  s_path/=partition;
  val_s_path->set(s_path);
  val_z_path->set(-(1./lambda)*std::log(partition));
    
  // derivatives
  for (pairiter it = neighpair.begin(); it != neighpair.end(); ++it) {
    ii = (*it).first;
    tmp=lambda*expdists[ii]*(s_path-(ii+1))/partition;
    s_path_ders[ii] = tmp;
    z_path_ders[ii] = expdists[ii]/partition;
  }
  for(int i = 0; i < coefficients.size(); ++i){
    s_der = 0.;
    z_der = 0.;
    for (pairiter it = neighpair.begin(); it != neighpair.end(); ++it) {
      ii = (*it).first;
      s_der += s_path_ders[ii] * numerators[ii][i];
      z_der += z_path_ders[ii] * numerators[ii][i];
    }
    setDerivative(val_s_path,i,s_der);
    setDerivative(val_z_path,i,z_der);
  }    
//  log.printf("CALCULATION DONE! \n");
}
    
///
/// this function updates the needed argument list
///
void FuncPathGeneral::prepare(){
  // neighbor list: rank and activate the chain for the next step 

  // neighbor list: if neigh_size<0 never sort and keep the full vector
  // neighbor list: if neigh_size>0  
  //                if the size is full -> sort the vector and decide the dependencies for next step 
  //                if the size is not full -> check if next step will need the full dependency otherwise keep this dependencies 

  if (neigh_size>0) {
    if (neighpair.size()==path_cv_values.size()) { // I just did the complete round: need to sort, shorten and give it a go
      // sort the values  
      sort(neighpair.begin(),neighpair.end(),pairordering());
      // resize the effective list
      neighpair.resize(neigh_size);
      log.printf("  NEIGH LIST NOW INCLUDE INDEXES: ");
      for(int i=0;i<neigh_size;++i)log.printf(" %i ",neighpair[i].first);log.printf(" \n");
    } else {
      if (int(getStep()) % int(neigh_stride/getTimeStep()) == 0) {
        log.printf(" Time %f : recalculating full neighlist \n",getStep()*getTimeStep());
        neighpair.resize(path_cv_values.size());  
        for(unsigned i=0;i<path_cv_values.size();++i)neighpair[i].first=i; 
      }
    } 
  }
    
  requestArguments(allArguments);
}

}
}
