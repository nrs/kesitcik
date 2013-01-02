

#include "pr1.h"
#include "prgl.h"
#include "section.h"
#include "glutMaster.h"
#include "glutWindow.h"
#include <fstream>
GlutMaster  * glutMaster;
MeshWindow  * firstWindow = 0;

using namespace std;
void read_input (char *path);


// a concrete class that allows cracking: \eps_ultimate>\eps_cracking
class concrete2: public material{
public:
  Real sig(Real eps);
  concrete2(Real a, Real b, Real c, Real d, Real e, Real f);
  concrete2();
private:
  void init(Real f_ck, Real eps_0, Real eps_cr1, Real eps_cr2, 
            Real eps_cu1, Real eps_cu2);
};


class trilinear: public material{
public:
  Real sig(Real eps);
  trilinear(Real a, Real b, Real c, Real d, Real e);
  trilinear();
private:
  void init(Real E, Real E_sh, Real eps_y, Real eps_sh, Real eps_u);
};


struct reinf
{
  Real x[2];
  Real d;
  reinf(Real a, Real b, Real dia){x[0] = a; x[1] = b; d=dia;}
};

void tsect_rein(supermesh &s, Real width, Real height, Real a, Real b, 
                vector<reinf> &r);

void gl_window_supermesh(supermesh &s);

vector<reinf> myrein;
Real flange_width = 1200;
Real beam_width = 300;
Real beam_height = 650;
Real flange_thickness = 150;
Real concrete_str = 20;
Real steel_eps_sh = 0.01;
Real steel_sig_y = 420;
Real steel_sig_ult = 600;
Real steel_eps_ult = 0.15;
Real steel_young = 2e5;
//Real clear_cover = 70;

bool reinforcement_defined = false;
bool input_file = false;


int main(int argc, char** argv) {

//  cout << argv[1] << endl;
  if (!argv[1]) {
    cout << "No input file specified. Reverting to defaults.\n";
  }else{
    read_input(argv[1]);
  }
// reinforcement diameter
  Real diameter = 28;
  // create the reinforcement array
  if (!reinforcement_defined){
    myrein.push_back(reinf(-110,50,diameter));
    myrein.push_back(reinf(110 ,50,diameter));
    myrein.push_back(reinf(0   ,50,diameter));
    myrein.push_back(reinf(110 ,113,diameter));
    myrein.push_back(reinf(-110,113,diameter));
  }

  supermesh s;
  // create the t-section
  tsect_rein(s,beam_width,beam_height,flange_thickness,flange_width,myrein);

  // readjust the center of gravity as in the question
  s.cg = node(0,270);
  char *basename;
  char dummy[] = "lol";

  if (input_file){
    basename=argv[1];
  }else{
    basename=dummy;
  }

  m_vs_phi(s,1000,0., basename,true);
//  interaction(s,500, basename,true);
//  sample_materials(s,0,0.16,1000,basename);
//  sample_materials(s,-0.0037,0,1000,basename);

//  gl_window_supermesh(s);

  return 0;

}


void reinforcement_layer (Real width, Real y,  Real dia, int amount){
  for (unsigned int i = 0; i < amount; i++){
    myrein.push_back(reinf(-1*width/2+i*width/(amount-1),y,dia));

  }

}

void read_input (char *path){
  cout << "Opening file: " << path << endl;

  std::ifstream infile(path);

  if (!infile) {cout<< "Error reading input file. Reverting to defaults.\n";return;}
  input_file=true;
  std::string line;
  while (std::getline(infile, line))
  {
    std::istringstream iss(line); string tok;
    Real a,b,c,d,e; string s;
    int f;
    while (getline( iss, tok, ' ' ))
    {
      if (tok == "rein"){
        if (!(iss >> a >> b >> c)){
          cout << "I/O Error\n"; exit(0);  
        } 
        reinforcement_defined = true;
        myrein.push_back(reinf(a,b,c)); break;
        cout << a <<" "<< b << " " << c << endl;
      } else if (tok == "rein_layer"){
        if (!(iss >> a >> b >> c >> f)){
          cout << "I/O Error\n"; exit(0);  
        } 
        reinforcement_defined = true;
        reinforcement_layer(a,b,c,f); break;
      } else if (tok == "beam_width"){
        if (!(iss >> a)){cout << "I/O Error\n"; exit(0);} 
        beam_width = a; break;
      } else if (tok == "flange_width"){
        if (!(iss >> a)){cout << "I/O Error\n"; exit(0);} 
        flange_width = a; break;
      } else if (tok == "beam_height"){
        if (!(iss >> a)){cout << "I/O Error\n"; exit(0);} 
        beam_height = a; break;
      } else if (tok == "flange_thickness"){
        if (!(iss >> a)){cout << "I/O Error\n"; exit(0);} 
        flange_thickness = a; break;
      } else if (tok == "steel_eps_sh"){
        if (!(iss >> a)){cout << "I/O Error\n"; exit(0);} 
        steel_eps_sh = a; break;
      } else if (tok == "steel_sig_y"){
        if (!(iss >> a)){cout << "I/O Error\n"; exit(0);} 
        steel_sig_y = a; break;
      } else if (tok == "steel_sig_ult"){
        if (!(iss >> a)){cout << "I/O Error\n"; exit(0);} 
        steel_sig_ult = a; break;
      } else if (tok == "steel_eps_ult"){
        if (!(iss >> a)){cout << "I/O Error\n"; exit(0);} 
        steel_eps_ult = a; break;
      } else if (tok == "conrete_str"){
        if (!(iss >> a)){cout << "I/O Error\n"; exit(0);} 
        concrete_str = a; break;
      } // else if (tok == "clear_cover"){
      //   if (!(iss >> a)){cout << "I/O Error\n"; exit(0);} 
      //   clear_cover = a; break;
      // } 


      
//      cout << "word " << tok << '\n';
    }
    // if (!(iss >>  s >> a >> b >> c)){
    //   cout << "I/O Error\n"; exit(0);  
    // } // error
    // cout << c << b << endl;
    // // process pair (a,b)
  }

  infile.close();

}


void gl_window_supermesh(supermesh &s){
  glutMaster   = new GlutMaster();    
  firstWindow  = new MeshWindow(glutMaster, 700, 700, 0, 0, "tsect"); 
  firstWindow->glsupermesh = &s;
  firstWindow->init();
  glutMaster->CallGlutMainLoop();
}

void tsect_rein(supermesh &s, Real width, Real height, Real a, Real b, 
                vector<reinf> &r)
{
  const int nsides=12;

  Real side_clearance = width/3;

  s.meshes = list<mesh> (2);

  rectangular_section(s.meshes.front(),-1*width/2,width/2,height-a,0);
  mesh dummy;
  rectangular_section(dummy,-1*b/2,b/2,height,height-a);
  s.meshes.front().add(dummy);

  for (vector<reinf>::iterator it=r.begin(); it!=r.end(); it++){
    mesh lol;
    circular_section(lol,it->d/2,it->x,nsides);
    s.meshes.back().add(lol);
  }


  s.meshes.front().subtract(s.meshes.back());

// Triangulation  
  s.meshes.front().triangulate_mesh(800);//50
  s.meshes.back().triangulate_mesh(20);//5

// Materials
  s.meshes.front().mat = //new concrete2;
    // new concrete2(25*0.85,
    //               2e-3,
    //               -0.0035, 
    //               0.00015,
    //               -0.020,
    //               0.080);
    new concrete2(concrete_str,2e-3,-0.0035, 0.0,-0.1,0.1);
  s.meshes.back().mat = new trilinear;

// Print info
  s.meshes.front().print_info(); cout<<endl;
  s.meshes.back().print_info();

// Color
  s.meshes.front().color_mesh(0,0,4);
  s.meshes.back().color_mesh(0,0,1);
  s.calc_cg();


}


///////////////////////
// Concrete2 Follows //
///////////////////////

Real concrete2::sig(Real eps){
  if (eps < par[2] || eps > par[3]){
    return 0.;
  }
  Real result;
  if (eps < 0 && eps >= par[1]){
    result = -1*par[0] * ((2*eps/par[1]) - (eps/par[1])*(eps/par[1])) ;
  }else if (eps <par[1] && eps>= par[2]) {
    result = -1*par[0]  + .15 * par[0] * (eps-par[1])/(par[2]-par[1]);
  }else{
    result = 0.;//eps*par[6];
  }
  return result;
}

// If e>0 and e<=b Then
// Return a*(2*e/b-(e/b)^2)
// ElseIf e>b and e<c Then
// Return a-.15*a*(e-b)/(c-b)


concrete2::concrete2(Real f_ck, Real eps_0, Real eps_cr1, Real eps_cr2, 
                     Real eps_cu1, Real eps_cu2)
{
  init(f_ck, eps_0, eps_cr1, eps_cr2, eps_cu1, eps_cu2);

}

concrete2::concrete2(){
  init(concrete_str,2e-3,-0.0035, 0.00015,-0.0035,0.00015);
}

void concrete2::init(Real f_ck, Real eps_0, Real eps_cr1, Real eps_cr2, 
                     Real eps_cu1, Real eps_cu2){
  par = std::vector<Real> (7);
  par[0]=(Real) fabs(f_ck);
  par[1]=(Real) -1*fabs(eps_0);
  par[2]=(Real) -1*fabs(eps_cr1); 
  par[3]=(Real)    fabs(eps_cr2); 
  par[4] = (Real) -1*fabs(eps_cu1); 
  par[5] = (Real)    fabs(eps_cu2); 
  epsult[0] = par[4];
  epsult[1] = par[5];
  description="Concrete 2";
  epsulttol[0] = epsult[0]-TOLERANCE; 
  epsulttol[1] = epsult[1]+TOLERANCE;
  par[6]=0.35*sqrt(f_ck)/eps_cr2; // tension slope
}

///////////////////////
// Steel2 Follows    //
///////////////////////


Real trilinear::sig(Real eps){
  if (eps < epsulttol[0] || eps > epsulttol[1])
    return 0.;
  Real result;
  Real epsp = fabs(eps);
  Real pos = eps > 0 ? 1. : -1.;

  if (epsp < par[2]){
    result = pos*epsp*par[0];

  }else if (epsp > par[2] && epsp < par[3]){
    result = pos*par[2]*par[0];
        
  }else if (epsp > par[3] && epsp < par[4]){
    result = pos*(par[1]*(epsp-par[3]) + par[2]*par[0]);
  }else{ result = 0;}
  return result;
}


trilinear::trilinear(Real E, Real E_sh, Real eps_y, Real eps_sh, Real eps_u){
  init(E, E_sh, eps_y, eps_sh, eps_u);
}

trilinear::trilinear(){
  init(steel_young, (steel_sig_ult-steel_sig_y)/(steel_eps_ult-steel_eps_sh), steel_sig_y/steel_young, steel_eps_sh, steel_eps_ult);
}

void trilinear::init(Real E, Real E_sh, Real eps_y, Real eps_sh, Real eps_u){
  par = std::vector<Real> (5);
  par[0]=(Real) E;
  par[1]=(Real) E_sh;
//    par[2]=(Real) eps_cu;
  par[2] = (Real) fabs(eps_y);
  par[3] = (Real) fabs(eps_sh); 
  par[4] = (Real) fabs(eps_u); 

  epsult[0] = (Real) -1*fabs(eps_u); 
  epsult[1] = (Real) fabs(eps_u); 
  description = "Trilinear steel";
  epsulttol[0] = epsult[0]-TOLERANCE; 
  epsulttol[1] = epsult[1]+TOLERANCE;
}
