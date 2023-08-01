/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SUTTG, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <stdexcept>
#include <fstream>
#include <string>
#include <queue>
#include <limits>
#include <sstream>

#include "XMLExpatParser.h"

#include "InputParser.H"

InputParser::InputParser( const std::string XMLinput,
			  const int verbose )
  : Xbounds(2,0.),Ybounds(2,0.),debug(verbose)
{
  //
  // parse the XML and save the tree
  //
  RefPtr<XMLParser> parser = new XMLExpatParser;
  xml_tree = parser->parse(XMLinput);

  //
  // get output verbosity
  //
  verbosity = xml_tree->get_child("RunSettings")->attributes()->get_int("verbosity",0);

  //
  // parse the desired tabulation 
  //
  RefPtr<XMLElement> all_tabulations = xml_tree->get_child("Tabulation");

  //First get defaults
  atomic_weight= all_tabulations->attributes()->get_double("A");
  atomic_number= all_tabulations->attributes()->get_double("Z");
  rho_ref = all_tabulations->attributes()->get_double("RRef");
  t_ref   = all_tabulations->attributes()->get_double("TRef");

  //
  // find the desired tabulation 
  //
  const std::string runtab = xml_tree->get_child("RunSettings")->attributes()->get("tabulation");
  if (verbosity > 0) std::cerr << "Tabulating using " << runtab << " rules" << std::endl;
  RefPtr<XMLElement> etab = all_tabulations->get_child(runtab);
  parseTabulation(etab);

  //
  // find the desired model
  //
  runmodel = xml_tree->get_child("RunSettings")->attributes()->get("model");
  if (verbosity > 0) std::cerr << "Tabulating using " << runmodel << " model" << std::endl;

  //
  // parse the model
  //
  RefPtr<XMLElement> eosmodel = xml_tree->get_child("EOSModel");
  parseEOSModel(runmodel,eosmodel);

}

void InputParser::parseTabulation( RefPtr<XMLElement> & etab )
{
  //
  // get the basics
  //
  tabtype = etab->attributes()->get("type");
  outputfilebase = etab->attributes()->get("basename");
  meshvars = etab->attributes()->get("meshvars","RX");
  logvars = etab->attributes()->get_int("logvars",1);
  numthreads = etab->attributes()->get_int("numthreads",0);
  nodealgorithm = etab->attributes()->get("algorithm","optimize");

  //
  // get tolerances
  //
  tolerance = etab->attributes()->get_double("tolerance");
  ftolerance = etab->attributes()->get_double("ftolerance",1.e-2);

  //
  // get sample density
  //
  bsamples = etab->attributes()->get_int("boundarySamples",2);
  rsamples = etab->attributes()->get_int("regionSamples",2);

  //
  // get number of nodes to add each time
  //
  addnodes = etab->attributes()->get_int("addNodes",1);

  //
  // get parameters for optimization algorithm
  //
  steperrorlevel = etab->attributes()->get_double("steperrorlevel",1.0);
  stepmultiplier = etab->attributes()->get_double("stepmultiplier",0.1);

  //
  // get the bounds -- density and temperature for RX mesh
  //                   pressure and temperature for PX mesh
  //
  if (meshvars[0] == 'R') {
    RefPtr<XMLElement> rb = etab->get_child("RBounds");
    Xbounds[0] = rb->attributes()->get_double("lower");
    Xbounds[1] = rb->attributes()->get_double("upper");
  }
  else if (meshvars[0] == 'P') {
    RefPtr<XMLElement> rb = etab->get_child("PBounds");
    Xbounds[0] = rb->attributes()->get_double("lower");
    Xbounds[1] = rb->attributes()->get_double("upper");
  }    
  else throw std::runtime_error("InputParser::parseTabulation: invalid mesh variables "+meshvars);

  RefPtr<XMLElement> tb = etab->get_child("TBounds");
  Ybounds[0] = tb->attributes()->get_double("lower");
  Ybounds[1] = tb->attributes()->get_double("upper");

}

double InputParser::getSharedParam( RefPtr<XMLElement> & mi,
				    const std::string name,
				    const std::string param_name,
				    const int param_loc,
				    const double param_scale )
{
  //
  // search for the shared parameter
  //
  int index = 0;
  RefPtr<XMLElement> spe;
  for (int i=0;i<mi->count_children();i++) {
    spe = mi->get_child(i);
    if (spe->label() == "DParam" && spe->attributes()->get("name") == name) break;
    index = i+1;
  }
  if (index >= mi->count_children()) {
    throw std::runtime_error("InputParser::getSharedParam: failed to find shared parameter named "+name);
  }

  //
  // attribute reference for the shared parameter
  //
  const RefPtr<XMLAttributeList> & pattr = spe->attributes();

  //
  // setup shared parameter data if it does not exist
  //
  if (sp.find(name) == sp.end()) {
    sp[name] = SharedParamData(pattr->get_double("value"));
  }

  //
  // add to the shared parameter list of copy locations
  //
  sp[name].names.push_back(param_name);
  sp[name].plocs.push_back(param_loc);
  sp[name].scales.push_back(param_scale);

  //
  // return the initial value of the shared parameter scaled appropriately
  //
  return sp[name].value*sp[name].scales.back();

}

void InputParser::parseEOSModel( const std::string mname,
				 RefPtr<XMLElement> & mi )
{
  //
  // check that all model names are unique and that only double params are shared
  //
  std::map<std::string,int> element_strings;
  for (int i=0;i<mi->count_children();i++) element_strings[mi->get_child(i)->label()]++;
  for (std::map<std::string,int>::iterator i=element_strings.begin();i!=element_strings.end();i++) {
    if (i->first == "SParam") throw std::runtime_error("InputParser::parseEOSModel -- Found shared string parameter");
    else if (i->first == "IParam") throw std::runtime_error("InputParser::parseEOSModel -- Found shared string parameter");
    else if (i->first != "DParam" && i->second > 1) throw std::runtime_error("InputParser::parseEOSModel -- Found multiple EOS models named: "+i->first);
  }

  //
  // Read all the needed models. This works by queuing up each needed
  // model name and parsing that entry in the element. Shared
  // parameters are loaded on demand from the requesting model into a
  // parameter holder with empty model string.
  //
  std::queue<std::string> mq;
  
  //
  // start with the top level model
  //
  mq.push(mname);

  while(!mq.empty()) {
    //
    // get model element from queued name
    //
    RefPtr<XMLElement> curmod = mi->get_child(mq.front());
    mq.pop();

    //
    // current model name
    //
    const std::string mlabel = curmod->label();
    
    if (debug > 0) std::cerr << "current model " << mlabel << std::endl;

    //
    // skip the model if it is already initialized
    //
    if (p.find(mlabel) != p.end()) continue;

    //
    // get the model type
    //
    const std::string mstring = curmod->attributes()->get("type");

    if (debug > 0) std::cerr << mlabel << " model type " << mstring << std::endl;

    //
    // get the number of parameters via the factory
    //
    std::vector<int> pnum(3,0);
    eosfactory.getParamCounts(mstring,pnum);

    //
    // initialize the parameters class for this model
    //
    p[mlabel] = new EOSParam( pnum[0], pnum[1], pnum[2] );

    //
    // set the model name
    //
    p[mlabel]->model_name = mstring;

    //
    // load the parameters for the model
    //
    for (int i=0;i<curmod->count_children();i++) {
      //
      // loop over each of the defined parameters
      //
      RefPtr<XMLElement> param = curmod->get_child(i);
    
      //
      // get the type of parameter
      //
      const std::string plabel = param->label();

      //
      // attributes for the parameter
      //
      const RefPtr<XMLAttributeList> & pattr = param->attributes();

      if (debug > 0) std::cerr << "child " << i << " type " << plabel << " name " << pattr->get("name","") << " value " << pattr->get("value","") << std::endl;

      //
      // Handle the case where an order is specified instead of a name
      // (valid only for string parameters).
      //
      const std::string pname = pattr->get("name","");
      std::size_t ploc;
      if (pname == "") {
        if (plabel != "SParam") throw std::runtime_error("InputParser::parseEOSModel -- missing name for "+plabel+" in model "+mlabel);
        ploc = pattr->get_int("order");
        
        //
        // expand string parameter vector if needed
        //
        if (ploc >= p[mlabel]->params_string.size()) p[mlabel]->params_string.resize(ploc+1);
      }
      else {
        ploc = eosfactory.getParamLocation(mstring,plabel,pname);
      }

      if (debug > 0) std::cerr << mlabel << " set " << plabel << " parameter " << ploc << " to ";
      if (plabel == "DParam") {
        //
        // handle shared parameters by catching if default max double is returned
        //
        p[mlabel]->params_double[ploc] = pattr->get_double("value",std::numeric_limits<double>::max());
        
        if (p[mlabel]->params_double[ploc] == std::numeric_limits<double>::max()) {
          p[mlabel]->params_double[ploc] = getSharedParam(mi,pattr->get("shared"),
                                                          mlabel,ploc,
							  pattr->get_double("scale",1.0));
        }

        if (debug > 0) std::cerr << p[mlabel]->params_double[ploc] << std::endl;
      }
      if (plabel == "IParam") {
        p[mlabel]->params_int[ploc]    = pattr->get_int("value");
        if (debug > 0) std::cerr << p[mlabel]->params_int[ploc] << std::endl;
      }
      if (plabel == "SParam") {
        p[mlabel]->params_string[ploc] = pattr->get("value");
        if (debug > 0) std::cerr << p[mlabel]->params_string[ploc] << std::endl;
      }

    }

    //
    // Check the string parameters. Ensure none are empty strings and
    // add string parameters to the queue as new model names
    //
    for (std::size_t i=0;i<p[mlabel]->params_string.size();i++) {
      if (p[mlabel]->params_string[i] == "") {
        std::ostringstream oss;
        oss << "InputParser::parserEOSModel -- Found empty string parameter in location " << i << " in model " << mlabel;
        throw std::runtime_error(oss.str());
      }
      else {
        mq.push(p[mlabel]->params_string[i]);
      }
    }
  }

}

void InputParser::copySharedParam( const std::string pname,
                                   const double val,
                                   Parameters & params )
{
  //
  // find shared parameter
  //
  std::map< std::string, SharedParamData >::iterator isp = sp.find(pname);
  if (isp == sp.end()) {
    throw std::runtime_error("InputParser::copySharedParam -- parameter name "+pname+" not found in shared parameters");
  }

  //
  // get reference to shared parameter data
  //
  const SharedParamData & spd = isp->second;

  //
  // copy shared parameter values to all desired locations
  //
  for (std::size_t i=0;i<spd.names.size();i++) {
    params[spd.names[i]]->params_double[spd.plocs[i]] = val*spd.scales[i];
  }

}
