///////////////////////////////////////////////////
//
//   Eddy Current Separator
//
//   This program is based on the following
//   libraries:
//   - ChronoEngine
//   - Irrlicht
// 
//  
// ------------------------------------------------ 
///////////////////////////////////////////////////
 
    
#include "physics/CHapidll.h" 
#include "physics/CHsystem.h"
#include "physics/CHconveyor.h"
#include "physics/CHbodyAuxRef.h"
#include "core/ChFileutils.h"
#include "irrlicht_interface/CHbodySceneNode.h"
#include "irrlicht_interface/CHbodySceneNodeTools.h" 
#include "irrlicht_interface/CHirrApp.h"
#include "core/ChRealtimeStep.h"
#include "core/ChMath.h"
#include "core/ChDistribution.h"
#include <irrlicht.h>
#include <fstream>
#include "unit_PYTHON/ChPython.h"
#include "unit_POSTPROCESS/ChPovRay.h"

// Use the namespace of Chrono

using namespace chrono;
using namespace postprocess;

// Use the main namespaces of Irrlicht
using namespace irr;
         
using namespace core;
using namespace scene; 
using namespace video;
using namespace io; 
using namespace gui;
//
using namespace std;

// Static values valid through the entire program (bad
// programming practice, but enough for quick tests)

double STATIC_flow = 100; 
double STATIC_speed = -1; //[m/s]

const double mu0 = 0.0000012566; //vacuum permability [Tm/A]  ****** Edo
const double drumspeed = 261; //[rad/s]
const double numberofpoles = 9;
const double intensity = 0.32; 
const double drumdiameter = 0.233; //0.30;
const double eta = 0.0000181; // Air drag coefficent [N*s/m^2]
double particles_dt;
double debris_number = 0;
double max_numb_particles = 100;


// conveyor constant
const double conveyor_length=0.830; //edo from CAD 12/06 //conveyor_length=0.800;//***ALEX, from CAD
const double conveyor_width=0.449; //***ALEX, from CAD, was 0.6
const double conv_thick = 0.05; // non troppo sottile, altrimenti non va la collision detection! Non importa se compenetra il cilindro.

const double ro=1.225;  //fluid density (air) [Kg/m^3]

// fence constant
const double fence_width = 0.02;
const double fence_height=0.2;
// bin constant 
const double y_posbin=-0.815;
const double binbase=0.02; // base thickness of the bin
const double bin_length=1; // accorciato da 3 a 1, ida
const double bin_width=1.5;
const double bin_height=0.2;
const double n = 2; // number of bins (values from 2 to 4), a regime devo importare il cad dei 3 contenitori, ida
// splitter constant
const double x_splitter1=0;
const double x_splitter2=0;
const double x_splitter3=0;
const double splitter_width=0.01;
const double splitter_height=0.4;
// hopper constant
const double xnozzlesize = 0.1;//0.2;
const double znozzlesize = 0.182; //**from CAD, nozzle width
const double ynozzlesize = 0.1;//0.5;
const double ynozzle = 0.01;
const double xnozzle = -conveyor_length/2+xnozzlesize/2+fence_width; //portato avanti****ida

const double densityMetal =  8400; // sn-pb //8900;//rame//1820 vecchia densità;
const double densityPlastic = 946;// polipropilene //900 vecchia densità;


// Coordinate systems with position and rotation of important items in the 
// simulator. These are initializad with constant values, but if loading the
// SolidWorks model, they will be changed accordingly to what is found in the CAD 
// file (see later, where the SolidWorks model is parsed). 

ChCoordsys<> conveyor_csys( ChVector<>(0, 0-conv_thick, 0) ) ; // default position
//ChCoordsys<> drum_csys    ( ChVector<>(conveyor_length/2, -(drumdiameter*0.5)-conv_thick/2,0) );  // default position
//ChCoordsys<> drum_csys    ( ChVector<>(conveyor_length, -(drumdiameter*0.5),0) );  // default position  ***Edo
ChCoordsys<> drum_csys    ( ChVector<>(0,0,0) );  // default position  ***Edo
ChCoordsys<> nozzle_csys  ( ChVector<>(xnozzle, ynozzle, 0) ); // default position
ChCoordsys<> Splitter1_csys  ( ChVector<>(conveyor_length/2+0.2, -(drumdiameter*0.5)-conv_thick/2,0) );  // default position
ChCoordsys<> Splitter2_csys  ( ChVector<>(conveyor_length/2+0.4, -(drumdiameter*0.5)-conv_thick/2,0) );  // default position
ChCoordsys<> Spazzola_csys  ( ChVector<>(conveyor_length/2-0.10, -(drumdiameter*0.5)-conv_thick/2,0) );  // default position


// set as true for saving log files each n frames
bool save_dataset = true;
bool save_irrlicht_screenshots = false;
bool save_POV_screenshots = false;
int saveEachNframes = 4;

bool irr_cast_shadows = true;

int totframes = 0;
	
bool init_particle_speed = false;

double particle_magnification = 1; // for larger visualization of particle ***Edo



// If one needs to add special data to a particle/object, the proper way to
// do this is by inheriting a custom class from the base ChAsset. 
// An asset of this class can then be created and attached to the
// particle (see below).  ***ALEX

class ElectricParticleProperty : public ChAsset
{ 
public:

	// data for this type of asset 
	ChVector<> Cdim;
	double conductivity;
	double birthdate;
	double chargeM;		//***ida + ale (coulomb, for plastic)
	double chargeP;
	ChVector<> InducedForce;
	ChVector<> InducedTorque;
	ChVector<> LiftForce;
	ChVector<> DragForce;
	enum 
	{
		e_fraction_box,
		e_fraction_cylinder,
		e_fraction_sphere,
		e_fraction_others
	} fraction;
	enum 
	{
		e_mat_plastic,
		e_mat_metal,
		e_mat_other
	} material_type;

	// default constructor with initialization
	ElectricParticleProperty()
	{
		Cdim = ChVector<>(1,1,1);
		conductivity = 0;
		birthdate = 0;
		fraction = ElectricParticleProperty::e_fraction_others;
		material_type = ElectricParticleProperty::e_mat_other;
		chargeM = 0;
		chargeP = 0;
	}
	
};





// This can be added to store the trajectory on a per-particle basis.

class ParticleTrajectory : public ChAsset
{
public:
	std::list< ChVector<> > positions;
	std::list< ChVector<> > speeds;
	unsigned int max_points;

	ParticleTrajectory() 
	{
		max_points = 80;
	}
};







// Define a MyEventReceiver class which will be used to manage input
// from the GUI graphical user interface

class MyEventReceiver : public IEventReceiver
{
public:

	MyEventReceiver(ChIrrAppInterface* myapp)
			{
				// store pointer applicaiton
				application = myapp;

				// ..add a GUI slider to control particles flow
				scrollbar_flow = application->GetIGUIEnvironment()->addScrollBar(
								true, rect<s32>(560, 15, 700, 15+20), 0, 101);
				scrollbar_flow->setMax(100);
				//scrollbar_flow->setPos(25);
				scrollbar_flow->setPos(((int)((STATIC_flow/1000.0 )*25)));
				text_flow = application->GetIGUIEnvironment()->addStaticText(
							L"Flow [particles/s]", rect<s32>(710,15,800,15+20), false);

				// ..add GUI slider to control the speed
				/*scrollbar_speed = application->GetIGUIEnvironment()->addScrollBar(
								true, rect<s32>(560, 40, 700, 40+20), 0, 102);
				scrollbar_speed->setMax(100); 
				//scrollbar_speed->setPos(72);
			    scrollbar_speed->setPos(((int)((STATIC_speed/2.0 )*67)));
				text_speed = application->GetIGUIEnvironment()->addStaticText(
								L"Conv.vel. [m/s]:", rect<s32>(710,40,800,40+20), false);*/ //deactivated for the moment!

				// ..add GUI checkmark to enable plotting forces
				checkbox_plotforces = application->GetIGUIEnvironment()->addCheckBox(false,core::rect<s32>(560,65, 560+150,65+20),
								0, 105, L"Plot applied CES forces");

				// ..add GUI checkmark to enable plotting forces
				checkbox_plottrajectories = application->GetIGUIEnvironment()->addCheckBox(false,core::rect<s32>(560,90, 560+150,90+20),
								0, 106, L"Plot trajectories");
				

			}

	bool OnEvent(const SEvent& event)
			{

				// check if user moved the sliders with mouse..
				if (event.EventType == EET_GUI_EVENT)
				{
					s32 id = event.GUIEvent.Caller->getID();
					IGUIEnvironment* env = application->GetIGUIEnvironment();
					switch(event.GUIEvent.EventType)
					{
					case EGET_SCROLL_BAR_CHANGED:
							if (id == 101) // id of 'flow' slider..
							{
								s32 pos = ((IGUIScrollBar*)event.GUIEvent.Caller)->getPos();
								STATIC_flow = 1000* ((double)pos/25);
								char message[50]; sprintf(message,"Flow %d [particl/s]", (int)STATIC_flow);
								text_flow->setText(core::stringw(message).c_str());
							}
							if (id == 102) // id of 'speed' slider..
							{
								s32 pos = ((IGUIScrollBar*)event.GUIEvent.Caller)->getPos();
								STATIC_speed = 3.0* (((double)pos)/100);
								char message[50]; sprintf(message,"Conv.vel. %2.2f [m/s]", STATIC_speed);
								text_speed->setText(core::stringw(message).c_str());
							}
					break;
					}
				} 
				// check if user presses keys
				if (event.EventType == irr::EET_KEY_INPUT_EVENT && !event.KeyInput.PressedDown)
				{
					irr::core::matrix4 matrix;
					switch (event.KeyInput.Key)
					{
					case irr::KEY_F1:	// camera will point to drum
						application->GetSceneManager()->getActiveCamera()->setPosition( vector3dfCH( drum_csys.pos + ChVector<>(0,0.4,0.4) ) );
						application->GetSceneManager()->getActiveCamera()->setTarget  ( vector3dfCH( drum_csys.pos ) );
						application->GetSceneManager()->getActiveCamera()->setFOV(36*chrono::CH_C_DEG_TO_RAD);
						return true;
					case irr::KEY_F2:	// camera will point to nozzle
						application->GetSceneManager()->getActiveCamera()->setPosition( vector3dfCH( nozzle_csys.pos + ChVector<>(0.2,0.3,0.4) ) );
						application->GetSceneManager()->getActiveCamera()->setTarget  ( vector3dfCH( nozzle_csys.pos ) );
						application->GetSceneManager()->getActiveCamera()->setFOV(36*chrono::CH_C_DEG_TO_RAD);
						return true;
					case irr::KEY_F3:	// camera will point to bins
						application->GetSceneManager()->getActiveCamera()->setPosition( vector3dfCH( drum_csys.pos + ChVector<>(0.8, 0.2, 0.9) ) );
						application->GetSceneManager()->getActiveCamera()->setTarget  ( vector3dfCH( drum_csys.pos + ChVector<>(0.2f,-0.8f, 0.f  ) ) );
						application->GetSceneManager()->getActiveCamera()->setFOV(36*chrono::CH_C_DEG_TO_RAD);
						return true;
					case irr::KEY_F4:	// camera will point to drum in orthographic projection
						application->GetSceneManager()->getActiveCamera()->setPosition( vector3dfCH( drum_csys.pos + ChVector<>(0.0,0.0, 4) ) );
						application->GetSceneManager()->getActiveCamera()->setTarget  ( vector3dfCH( drum_csys.pos ) );	
						matrix.buildProjectionMatrixOrthoLH(0.4f, 0.3f, 0.3f, 100.f);
						application->GetSceneManager()->getActiveCamera()->setProjectionMatrix(matrix,true);
						return true;

					default:
						break;
					}
				}

				return false;
			}

private: 
	ChIrrAppInterface* application;

public:
	IGUIScrollBar*  scrollbar_flow;
	IGUIStaticText* text_flow;
	IGUIScrollBar*  scrollbar_speed;
	IGUIStaticText* text_speed;
	IGUICheckBox*	checkbox_plotforces;
	IGUICheckBox*	checkbox_plottrajectories;
};


// Shortcut for easy creation of a body that has both a box visualization and a box collision,
// with default values for mass inertia etc. 

ChSharedPtr<ChBody> create_box_collision_shape(ChVector<> pos, //  center of box
											   ChVector<> size, // size along box x y z
											   ChQuaternion<> rot = QUNIT) // optional rotation
{
	// Create a body
	ChSharedPtr<ChBody> mrigidBody(new ChBody);

	mrigidBody->SetPos(pos);
			
	// Define a collision shape 
	mrigidBody->GetCollisionModel()->ClearModel();
	mrigidBody->GetCollisionModel()->AddBox(0.5*size.x, 0.5*size.y, 0.5*size.z);
	mrigidBody->GetCollisionModel()->BuildModel();
	mrigidBody->SetCollide(true);

	// Attach a visualization shape asset. 
	ChSharedPtr<ChBoxShape> mbox(new ChBoxShape);
	mbox->GetBoxGeometry().SetLenghts(size);
	mrigidBody->AddAsset(mbox);

	return mrigidBody;
}


// Function that creates debris that fall on the conveyor belt, to be called at each dt

void create_debris(double dt, double particles_second, 
				   ChCoordsys<> mnozzle_csys,
				   ChSystem& mysystem, 
				   ChIrrApp* irr_application,
				   ChPovRay* mpov_exporter)
{

	double sph_fraction = 0.7;
	double box_fraction = 0.3;
	double cyl_fraction = 1-box_fraction-sph_fraction;

	//double sphrad = 0.6e-3 + (ChRandom()-0.5)*(0.6e-3); vecchia distribuzione
	double sphrad = 0.009;//3e-3; 
	double cylhei = 0.035;
	double cylrad = sphrad;
	double cylmass = CH_C_PI*pow(cylrad,2)*cylhei* 1.0;  // now with default 1.0 density
	double sphmass = (4./3.)*CH_C_PI*pow(sphrad,3)* 1.0; // now with default 1.0 density
	double sphinertia = 0.4*pow(sphrad,2)*sphmass;      // now with default 1.0 density
	double cylinertia = 0.0833*(pow(cylhei,2)+3*pow(sphrad,2))*cylmass;//0.0833*(pow(cylhei,2)+3*pow(sphrad,2))*cylmass;  // now with default 1.0 density
	double cylinertia2 = 0.5*pow(sphrad,2)*cylmass; //0.5*pow(sphrad,2)*cylmass; // now with default 1.0 density
	

	double exact_particles_dt = dt * particles_second;
	//double particles_dt= floor(exact_particles_dt);
	particles_dt= floor(exact_particles_dt);
	double remaind = exact_particles_dt - particles_dt;

    
	if (remaind > ChRandom()) particles_dt +=1;

    
	for (int i = 0; i < particles_dt; i++)
	{
		ChSharedPtr<ChBody> created_body;
		ChSharedPtr<ElectricParticleProperty> created_electrical_asset;

		double rand_shape_fract = ChRandom();

		
		ChVector<> rand_position = mnozzle_csys.TrasformLocalToParent( ChVector<>(-0.35*xnozzlesize+ChRandom()*xnozzlesize, 0, -0.5*znozzlesize+ChRandom()*znozzlesize) ) ;
     
		//
		// 1 ---- Create particle 
		// 

		if (rand_shape_fract < sph_fraction)
		{
			// Create a body
			ChSharedPtr<ChBody> mrigidBody(new ChBody);

			mrigidBody->SetPos(rand_position);
			mrigidBody->SetMass(sphmass);
			mrigidBody->SetInertiaXX(ChVector<>(sphinertia,sphinertia,sphinertia));
			mrigidBody->SetFriction(0.2f);
		    mrigidBody->SetRollingFriction(0.2f);
			mrigidBody->SetNoGyroTorque(true);
			
			

			// Define a collision shape 
			mrigidBody->GetCollisionModel()->ClearModel();
			mrigidBody->GetCollisionModel()->AddSphere(sphrad);
			//mrigidBody->GetCollisionModel()->AddSphere(sphrad, &ChVector<>(0.005,0,0)); // etc. optional
			mrigidBody->GetCollisionModel()->BuildModel();
			mrigidBody->SetCollide(true);

			// Attach a visualization shape asset. 
			ChSharedPtr<ChSphereShape> msphere(new ChSphereShape);
			msphere->GetSphereGeometry().rad = sphrad * particle_magnification; 
			msphere->GetSphereGeometry().center = ChVector<>(0,0,0);
			mrigidBody->AddAsset(msphere);
			/* etc. optional
			ChSharedPtr<ChSphereShape> msphere2(new ChSphereShape);
			msphere2->GetSphereGeometry().rad = sphrad*5; // test****
			msphere2->GetSphereGeometry().center = ChVector<>(0.005,0,0);
			mrigidBody->AddAsset(msphere2);
			*/

			// Attach a custom asset. Look how to create and add a custom asset to the object! ***ALEX
			ChSharedPtr<ElectricParticleProperty> electric_asset(new ElectricParticleProperty); 
			electric_asset->Cdim         = ChVector<>(2*sphrad,2*sphrad,2*sphrad);
			electric_asset->fraction	 = ElectricParticleProperty::e_fraction_sphere;
			electric_asset->birthdate	 = mysystem.GetChTime();
			mrigidBody->AddAsset(electric_asset);
			
			// Finally, do not forget to add the body to the system:
			mysystem.Add(mrigidBody);

			created_body = mrigidBody;
			created_electrical_asset = electric_asset; // for reference later
		}//**********************************************************************************************************************
		
		if ((rand_shape_fract > sph_fraction) && 
			(rand_shape_fract < box_fraction+sph_fraction))
		{
			double xscale = 1.3*(1-0.8*ChRandom()); // for oddly-shaped boxes..
			double yscale = 1.3*(1-0.8*ChRandom());
			double zscale = 1.3*(1-0.8*ChRandom());

			//	ChCollisionModel::SetDefaultSuggestedEnvelope(0.002);
			//	ChCollisionModel::SetDefaultSuggestedMargin  (0.002);

			// Create a body
			ChSharedPtr<ChBody> mrigidBody(new ChBody);

			mrigidBody->SetPos(rand_position);
			mrigidBody->SetMass(sphmass);
			mrigidBody->SetInertiaXX(ChVector<>(sphinertia,sphinertia,sphinertia));
			mrigidBody->SetFriction(0.4f);
			mrigidBody->SetImpactC(0.0f); 
			mrigidBody->SetNoGyroTorque(true);

			// Define a collision shape 
			mrigidBody->GetCollisionModel()->ClearModel();
			mrigidBody->GetCollisionModel()->AddBox(sphrad*2*xscale, sphrad*2*yscale, sphrad*2*yscale);
			mrigidBody->GetCollisionModel()->BuildModel();
			mrigidBody->SetCollide(true);

			// Attach a visualization shape asset. 
			ChSharedPtr<ChBoxShape> mbox(new ChBoxShape);
			mbox->GetBoxGeometry().SetLenghts(ChVector<>(sphrad*2*xscale, sphrad*2*yscale, sphrad*2*yscale) * particle_magnification);
			mrigidBody->AddAsset(mbox);

			// Attach a custom asset. Look how to create and add a custom asset to the object! ***ALEX
			ChSharedPtr<ElectricParticleProperty> electric_asset(new ElectricParticleProperty); 
			electric_asset->Cdim         = ChVector<>(2*sphrad,2*sphrad,2*sphrad);
			electric_asset->fraction	 = ElectricParticleProperty::e_fraction_box;
			electric_asset->birthdate	 = mysystem.GetChTime();
			mrigidBody->AddAsset(electric_asset);
			
			// Finally, do not forget to add the body to the system:
			mysystem.Add(mrigidBody);

			created_body = mrigidBody; // for reference later
			created_electrical_asset = electric_asset; // for reference later
		}

		if (rand_shape_fract > box_fraction+sph_fraction)
		{

			//	ChCollisionModel::SetDefaultSuggestedEnvelope(0.002);
			//	ChCollisionModel::SetDefaultSuggestedMargin  (0.002);

			// Create a body
			ChSharedPtr<ChBody> mrigidBody(new ChBody);

			mrigidBody->SetPos(rand_position);
			mrigidBody->SetMass(cylmass);
			mrigidBody->SetInertiaXX(ChVector<>(cylinertia,cylinertia2,cylinertia));
			mrigidBody->SetFriction(0.4f);
			mrigidBody->SetImpactC(0.0f);
			mrigidBody->SetNoGyroTorque(true);
			

			// Define a collision shape 
			mrigidBody->GetCollisionModel()->ClearModel();
			mrigidBody->GetCollisionModel()->AddCylinder(sphrad,sphrad,cylhei);
			mrigidBody->GetCollisionModel()->BuildModel();
			mrigidBody->SetCollide(true);

			// Attach a visualization shape asset. 
			ChSharedPtr<ChCylinderShape> mcyl(new ChCylinderShape);
			mcyl->GetCylinderGeometry().rad = sphrad;
			mcyl->GetCylinderGeometry().p1 = ChVector<>(0, cylhei/2,0) * particle_magnification;
			mcyl->GetCylinderGeometry().p2 = ChVector<>(0,-cylhei/2,0) * particle_magnification;
			mrigidBody->AddAsset(mcyl);

			// Attach a custom asset. Look how to create and add a custom asset to the object! ***ALEX
			ChSharedPtr<ElectricParticleProperty> electric_asset(new ElectricParticleProperty); 
			electric_asset->Cdim         = ChVector<>(sphrad*2,cylhei,sphrad*2);
			electric_asset->fraction	 = ElectricParticleProperty::e_fraction_cylinder;
			electric_asset->birthdate	 = mysystem.GetChTime();
			mrigidBody->AddAsset(electric_asset);
			
			// Finally, do not forget to add the body to the system:
			mysystem.Add(mrigidBody);

			created_body = mrigidBody; // for reference later
			created_electrical_asset = electric_asset; // for reference later
		}

		//
		// 2 ---- Adjust stuff that depends on the metal/plastic fraction
		//

		// Depending on a randomized fraction, set the material type as 'plastic' 
		// or 'metal', by setting the 'material_type' in the ElectricParticleProperty of the
		// created particle and the 'conductivity' (regardless if it was a cylinder 
		// or cube etc.) ***ALEX
		if (!created_electrical_asset.IsNull())
		{
			  
			double rand_mat = ChRandom();

			double plastic_fract = 0.5;
				
			if (rand_mat < plastic_fract)
			{
				created_electrical_asset->conductivity = 0;
                created_electrical_asset->material_type = ElectricParticleProperty::e_mat_plastic;

				// Attach a 'pink' texture to easily recognize plastic stuff in 3d view
				ChSharedPtr<ChTexture> mtexture(new ChTexture);
				mtexture->SetTextureFilename("../objects/pinkwhite.png");
				created_body->AddAsset(mtexture);

				// Multiply the default mass & intertia tensor by density (previously assumed =1)
				
				created_body->SetDensity(densityPlastic);
				created_body->SetMass( created_body->GetMass() * ::densityPlastic);
				created_body->SetInertiaXX( created_body->GetInertiaXX() * ::densityPlastic);
			} 
			if (rand_mat > plastic_fract)
			{
				created_electrical_asset->conductivity = 6428000;//6670000 conducibilità vecchia;
				created_electrical_asset->material_type = ElectricParticleProperty::e_mat_metal;

				// Attach a 'blue' texture to easily recognize metal stuff in 3d view
				ChSharedPtr<ChTexture> mtexture(new ChTexture);
				mtexture->SetTextureFilename("../objects/bluwhite.png");
				created_body->AddAsset(mtexture);

				// Multiply the default mass & intertia tensor by density (previously assumed =1)
		
				created_body->SetDensity(densityMetal);
				created_body->SetMass( created_body->GetMass() * ::densityMetal);
                created_body->SetInertiaXX( created_body->GetInertiaXX() * ::densityMetal);
			}


			if (created_electrical_asset->fraction == ElectricParticleProperty::e_fraction_sphere)
			{
			}

		}

		//
		// 3 ---- Set parameters that are common for whatever created particle, ex. limit speed threshold:
		//

		// This is an optional hack that largely affects the stability of the
		// simulation. 
		// In fact, if particles happen to spin too fast, the collision detection
		// starts to be very slow, and maybe also inaccurate; also, the time integration
		// could diverge. To get rid of this problem wihtout reducing too much the timestep, 
		// one can enable a limit on angular velocity and/or linear velocity. NOTE that 
		// this achieves greater stability at the cost of lower realism of the simulation, 
		// so it should not be abused. ***ALEX

		bool do_velocity_clamping = true;

		if (!created_body.IsNull() && do_velocity_clamping)
		{
			created_body->SetLimitSpeed(true);
			created_body->SetMaxSpeed(100);
			created_body->SetMaxWvel(250);
		}

		//
		// 4 ---- Irrlicht setup for enabling the 3d view of the particle
		//

		// If Irrlicht is used, setup also the visualization proxy.
		// This is necessary so that Irrlicht 3D view can show the visualization
		// assets that have been added to the created bodies. **ALEX
		if (!created_body.IsNull() && irr_application)
		{
			irr_application->AssetBind(created_body);
			irr_application->AssetUpdate(created_body);
		}

		//
		// 5  ---- If a POV exporter is used, one must add the particle to its database
		// 

		if (!created_body.IsNull() && mpov_exporter && save_POV_screenshots)
		{
			mpov_exporter->Add(created_body);
		}


		// 
		// 5 ---- it could be useful to attach an asset for storing the trajectory of particle
		//

		if (!created_body.IsNull())
		{
			ChSharedPtr<ParticleTrajectory> massettraj(new ParticleTrajectory);
			created_body->AddAsset(massettraj);
		}

		// 
		// 6 --- set initial speed to particles as soon as created, so they do not have transient slip with belt
		//

		if (init_particle_speed)
		{
			created_body->SetPos_dt(ChVector<>(0,0,0));
			//created_body->SetPos_dt(ChVector<>(STATIC_speed, 0,0));
		}

	}

} 


     
  
// Function that deletes old debris (to avoid infinite creation that fills memory)

void purge_debris(ChSystem& mysystem, double max_age = 5.0)
{
	for (unsigned int i=0; i<mysystem.Get_bodylist()->size(); i++)
	{
		ChBody* abody = (*mysystem.Get_bodylist())[i];

		bool to_delete = false;

		// Fetch the ElectricParticleProperty asset from the list of 
		// assets that have been attached to the object, and retrieve the
		// custom data that have been stored. ***ALEX
		for (unsigned int na= 0; na< abody->GetAssets().size(); na++)
		{
			ChSharedPtr<ChAsset> myasset = abody->GetAssetN(na);
			if (myasset.IsType<ElectricParticleProperty>())
			{
				ChSharedPtr<ElectricParticleProperty> electricproperties = myasset.DynamicCastTo<ElectricParticleProperty>();
				double particle_birthdate  = electricproperties->birthdate ;
				double particle_age = mysystem.GetChTime() - particle_birthdate;
				if (particle_age > max_age)
				{
					to_delete = true;
				}
			}
		}

		if (to_delete)
		{
			abody->AddRef();	// dirty trick to convert basic pointer to..
			ChSharedPtr<ChBody> mysharedbody(abody); // ..shared pointer

			mysystem.Remove(mysharedbody);
			// mysharedbody->RemoveRef(); //***NOT needed - previously needed cause always Add() to POV exporter..
			i--; // this because if deleted, the rest of the array is shifted back one position..
		}
	}

}

// ECS forces *****Edo

void apply_forces (	ChSystem* msystem,		// contains all bodies
						ChCoordsys<>& drum_csys, // pos and rotation of drum 
						double drumspeed,		 // speed of drum
						double numberofpoles,	 // number of couples of poles
						double intensity,		 // intensity of the magnetic field
						double drumdiameter,
						int totframes)		
{
	
		for (unsigned int i=0; i<msystem->Get_bodylist()->size(); i++)
	{
		ChBody* abody = (*msystem->Get_bodylist())[i];

		bool was_a_particle = false;
		ChSharedPtr<ElectricParticleProperty> electricproperties; // null by default

		// Fetch the ElectricParticleProperty asset from the list of 
		// assets that have been attached to the object, and retrieve the
		// custom data that have been stored. ***ALEX
		for (unsigned int na= 0; na< abody->GetAssets().size(); na++)
		{
			ChSharedPtr<ChAsset> myasset = abody->GetAssetN(na);
			if (myasset.IsType<ElectricParticleProperty>())
			{
				// OK! THIS WAS A PARTICLE! ***ALEX
				was_a_particle = true;	
				electricproperties = myasset.DynamicCastTo<ElectricParticleProperty>();
				//electricproperties = myasset;
			} 
		}

		// Do the computation of forces only on bodies that had 
		// the 'ElectricParticleProperty' attached.. **ALEX
		if(was_a_particle)
		{

			ChVector<> diam = electricproperties->Cdim; 
			double sigma =    electricproperties->conductivity;

			// Remember to reset 'user forces accumulators':
		abody->Empty_forces_accumulators();

		// initialize speed of air (steady, if outside fan stream): 
		ChVector<> abs_wind(0,0,0);
 
		// calculate the position of body COG with respect to the drum COG:
		
		ChVector<> mrelpos = drum_csys.TrasformParentToLocal(abody->GetPos());
		double distx=(mrelpos.x-0.246771402771896);
		double disty=(mrelpos.y+0.53927223197952);
		ChVector<> velocity=abody->GetPos_dt();
		double velocityx=velocity.x;
		double velocityy=velocity.y;
		double velocityz=velocity.z;
		ChVector <> rot_speed=abody->GetWvel_par();
		double rot_speedz=rot_speed.z; //bisogna tirare fuori la componente attorno all'asse z globale della velocità di rotazione

		double velocity_norm_sq=velocity.Length2();

		//ChQuaternion<> rot_velocity=abody->GetRot_dt;
		double distance = pow(distx*distx+disty*disty,0.5);//mrelpos.Length();
		double phi = atan2(disty,distx);
		double phi2 = atan2(-velocity.y,velocity.x);
		double B = intensity*pow((drumdiameter/2/distance),numberofpoles+1);
		double shape = abody->GetCollisionModel()->GetShapeType();
		ChMatrix33<> RotM = *abody->GetA();
		ChVector<> yG = RotM.Matr_x_Vect(VECT_Y);
		ChVector<> particleRspeed = abody->GetWvel_par();
		ChVector<> particleRspeed_loc = abody->GetWvel_loc();
		double rho = acos(yG.y/yG.Length());
		
		

		if (rho > CH_C_PI_2)
		{
			rho = CH_C_PI-rho;
		}
		double alpha = rho/CH_C_PI_2;

		// check velocità rotazione particelle
	
		double constR;
		double constI;
		//double constTorque;
		double Volume,d;

		
		if (shape == 0) //sphere
		{
			d=diam.x;
			double csi = mu0*d*d* sigma *(drumspeed*numberofpoles+particleRspeed.z);
			constR = 21*pow(csi,2)/(20*(1764+pow(csi,2)));
			constI = 21*42*pow(csi,2)/(20*(1764+pow(csi,2)));
			//constTorque = 1/40*pow(diam.x,2)*(4/3*CH_C_PI*pow(diam.x/2,3))* sigma*(drumspeed*numberofpoles+particleRspeed.z);
			Volume = (4./3.)*CH_C_PI*pow(diam.x/2,3);
			//constTorque = -pow(B,2)*4/3*CH_C_PI*pow(diam.x/2,3)*constI/mu0;


			//data_forces << csi << ",\t\t\t";
		}
		else if (shape == 2) // box (as a sphere)
		{
			d=diam.x;
			//double csi = mu0*pow(diam.x,2)* sigma*(drumspeed*numberofpoles+particleRspeed.z);
			//constR = 21*pow(csi,2)/20/(1764+pow(csi,2));
			//constI = 21*42*pow(csi,2)/(1764+pow(csi,2));
			double csi = mu0*d*d* sigma *(drumspeed*numberofpoles+particleRspeed.z);
			constR = 21*pow(csi,2)/(20*(1764+pow(csi,2)));
			constI = 21*42*pow(csi,2)/(20*(1764+pow(csi,2)));
			//constTorque = 1/40*pow(diam.x,2)*(4/3*CH_C_PI*pow(diam.x/2,3))* sigma*(drumspeed*numberofpoles+particleRspeed.z);
			Volume = pow(diam.x,3);
			//constTorque = -pow(B,2)*pow(diam.x,3)*constI/mu0;

			//data_forces << csi << ",\t\t\t";
		}
		else if (shape == 3) // cylinder
		{
			if (diam.y<diam.x) // disk
			{
				d=diam.x; //se togliamo alpha=0, cambiarlo...
				double csi_par = mu0*pow(diam.y,2)* sigma*(drumspeed*numberofpoles+particleRspeed.z);
				double csi_per = mu0*pow(diam.x,2)* sigma*(drumspeed*numberofpoles+particleRspeed.z);
				double den = 4*(256+pow(0.6*CH_C_PI*diam.y*csi_per/diam.x,2));
				constR = alpha*0.6*CH_C_PI*diam.y*pow(csi_per,2)/diam.x/den+(1-alpha)*pow(csi_par,2)/(144+pow(csi_par,2));
				constI = alpha*16*csi_per/den+(1-alpha)*12*csi_par/(144+pow(csi_par,2));
				//constTorque = (alpha*1/64*(diam.x,2)+(1-alpha)*1/12*(diam.y,2))* sigma*(drumspeed*numberofpoles+particleRspeed.z)* (CH_C_PI*pow(diam.x/2,2)*diam.y);
				Volume = CH_C_PI*pow(diam.x/2,2)*diam.y;
				//constTorque = -pow(B,2)*CH_C_PI*pow(diam.x/2,2)*diam.y*constI/mu0;
				//data_forces << csi_par << ",\t";
				//data_forces << csi_per << ",\t";
			}
			else // cylinder
			{	
				d=diam.x;
				double csi = mu0*pow(diam.x,2)* sigma*(drumspeed*numberofpoles+particleRspeed.z);
				constR = alpha*9*pow(csi,2)/(8*(576+pow(csi,2)))+(1-alpha)*3*pow(csi,2)/(2*(576+pow(csi,2)));
				constI = alpha*9*24*csi/(8*(576+pow(csi,2)))+(1-alpha)*3*24*csi/(2*(576+pow(csi,2)));
				//constTorque = (alpha*3/64+(1-alpha)*1/16)*pow(diam.x,2)* sigma*(drumspeed*numberofpoles+particleRspeed.z)* (CH_C_PI*pow(diam.x/2,2)*diam.y);
				Volume = CH_C_PI*pow(diam.x/2,2)*diam.y;
				//constTorque = -pow(B,2)*CH_C_PI*pow(diam.x/2,2)*diam.y*constI/mu0;

				//data_forces << csi << ",\t\t\t";
			}
		}
		double R0=(rot_speedz*2*d)/(2*pow(velocity_norm_sq,0.5));
		double CL=R0/(2.2*R0+0.7);
	    double CLdisk=1.4*R0/(R0+1);
		double CD=pow(0.1+2*pow(CL,2),0.5);
	    double CDdisk=pow(0.16+0.25*pow(CLdisk,2),0.5);
		
		if (shape==3 && diam.y<diam.x)
		{
			CL=CLdisk;
			CD=CDdisk;
		}

	

		ChVector<> DragForce = electricproperties->DragForce;
		electricproperties->DragForce.x = -CD*ro*velocity_norm_sq*CH_C_PI*diam.x*diam.y/2*cos(phi2);
		electricproperties->DragForce.y = CD*ro*velocity_norm_sq*CH_C_PI*diam.x*diam.y/2*sin(phi2);
		electricproperties->DragForce.z = 0;

		abody->Accumulate_force( DragForce, abody->GetPos(), false);

		ChVector<> LiftForce = electricproperties-> LiftForce;
		electricproperties->LiftForce.x = CL*ro*velocity_norm_sq*CH_C_PI*diam.x*diam.y/2*sin(phi2);
	    electricproperties->LiftForce.y = CL*ro*velocity_norm_sq*CH_C_PI*diam.x*diam.y/2*cos(phi2);
		electricproperties->LiftForce.z = 0;	
	
		abody->Accumulate_force(LiftForce, abody->GetPos(), false);

	    //if (distance>0.5)
		//{
		ChVector<> InducedForce = electricproperties->InducedForce;
		/*double Induced_Fr = ((numberofpoles+1)*pow(B,2)*Volume/mu0/distance)*constR;
		double Induced_Fphi = ((numberofpoles+1)*pow(B,2)*Volume/mu0/distance)*constI;
		double InducedF = sqrt(pow(Induced_Fr,2)+pow(Induced_Fr,2));*/
		electricproperties->InducedForce.x = ((numberofpoles+1)*pow(B,2)*Volume/mu0/distance)*(constR*cos(phi)+constI*sin(phi));
		electricproperties->InducedForce.y = ((numberofpoles+1)*pow(B,2)*Volume/mu0/distance)*(constR*sin(phi)-constI*cos(phi));
		electricproperties->InducedForce.z = 0;	
		abody->Accumulate_force(InducedForce, abody->GetPos(), false);
		//}

		ChVector<> InducedTorque = electricproperties->InducedTorque;
		electricproperties->InducedTorque.x = 0;
		electricproperties->InducedTorque.y = 0;
		//InducedTorque.z = -constTorque*pow(B,2);
		electricproperties->InducedTorque.z = (-pow(B,2)*Volume*constI)/mu0;
		abody->Accumulate_torque(InducedTorque, false);

		
		//coordinate del rotore. la y del rotore è la z delle coordinate del sistema
		ChVector<> pos = drum_csys.TrasformParentToLocal(abody->GetPos());
		ChVector<> acc_force=abody->Get_accumulated_force();
		ChVector<> acc_torque=abody->Get_accumulated_torque();
		
		ChVector<> iner=abody->GetInertiaXX();
		double ingxx = iner.x;
		double ingyy = iner.y;
		double ingzz = iner.z;

		double posx=pos.x;
		double posy=pos.y;
        
	    }
		
	}
}
 
// Function that defines the forces on the debris ****ida

/*void apply_forces (	ChSystem* msystem,		// contains all bodies
						ChCoordsys<>& drum_csys, // pos and rotation of drum 
						double drumspeed,		 // speed of drum
					    double drumdiameter,
						double h1,
						double h2,
						double L,
						double electrodediameter,
						double j,
						double alpha,
			            double U,
						double f,
					


						int totframes)		
{

	
	for (unsigned int i=0; i<msystem->Get_bodylist()->size(); i++)
	{
		ChBody* abody = (*msystem->Get_bodylist())[i];

		bool was_a_particle = false;
		ChSharedPtr<ElectricParticleProperty> electricproperties; // null by default

		// Fetch the ElectricParticleProperty asset from the list of 
		// assets that have been attached to the object, and retrieve the
		// custom data that have been stored. ***ALEX
		for (unsigned int na= 0; na< abody->GetAssets().size(); na++)
		{
			ChSharedPtr<ChAsset> myasset = abody->GetAssetN(na);
			if (myasset.IsType<ElectricParticleProperty>())
			{
				// OK! THIS WAS A PARTICLE! ***ALEX
				was_a_particle = true;		
				electricproperties = myasset.DynamicCastTo<ElectricParticleProperty>();
			} 
		}

		// Do the computation of forces only on bodies that had 
		// the 'ElectricParticleProperty' attached.. **ALEX
		if(was_a_particle)
		{

			ChVector<> diam = electricproperties->Cdim; 
			double sigma =    electricproperties->conductivity;

			// Remember to reset 'user forces accumulators':
			abody->Empty_forces_accumulators();

			// initialize speed of air (steady, if outside fan stream): 
			ChVector<> abs_wind(0,0,0);

			// calculate the position of body COG with respect to the drum COG:
			ChVector<> mrelpos = drum_csys.TrasformParentToLocal(abody->GetPos());
			double distx=mrelpos.x;
			double disty=mrelpos.y;
			ChVector<> velocity=abody->GetPos_dt();
			double velocityx=velocity.x;
			double velocityy=velocity.y;
			double velocityz=velocity.z;
			//ChVector <> rot_speed=abody->GetWvel_par();
			//double rot_speedz=rot_speed.z; //bisogna tirare fuori la componente attorno all'asse z globale della velocità di rotazione

			double velocity_norm_sq=velocity.Length2();

			//ChQuaternion<> rot_velocity=abody->GetRot_dt;
			
			// Polar coordinates of particles respect to the axis of the rotor, may be useful later **ALEX

			double distance = pow(distx*distx+disty*disty,0.5);
			double phi = atan2(disty,distx);
			double phi2 = atan2(-velocity.y,velocity.x);
			

		

			//
			// STOKES FORCES
			//


			double average_rad = 0.5* electricproperties->Cdim.Length(); // Approximate to sphere radius. Ida: this can be improved, by having Stokes forces for three Cdim x y z values maybe 

			ChVector<> StokesForce = electricproperties->StokesForce;
			
			electricproperties->StokesForce =0;
	
			abody->Accumulate_force(StokesForce, abody->GetPos(), false);


            

            //Calculating the analytical expressions of the electric field***ida

			double xuno=distx*cos(alpha)+ disty*sin(alpha);//analytical parameter****ida
            double yuno=disty*cos(alpha)- distx*sin(alpha);//analytical parameter****ida

			double Ex=(((j-h1+xuno)/(pow((j-h1+xuno),2)+pow(yuno,2))+((j+h1-xuno)/(pow((j+h1-xuno),2)+pow(yuno,2)))*f));//analytical expression of the electric field x direction***ida
			double Ey=((yuno/(pow((j-h1+xuno),2)+pow(yuno,2))-(yuno/(pow((j+h1-xuno),2)+pow(yuno,2)))*f));//analytical expression of the electric field y direction***ida
			double Ez=0;

			ChVector<> vE (Ex, Ey, Ez);
			double E = vE.Length();

          
			
		    //
			//===== METAL FORCES ==========
			//

			if (electricproperties->material_type == ElectricParticleProperty::e_mat_metal)
			{ 
				// charge the particle? (contact w. drum)
				if ((distx > 0) && (disty > 0))
				{
					if (electricproperties->chargeM == 0)
                    {
					electricproperties->chargeM = 0.666666666666667*pow(CH_C_PI,3)*epsilon*pow(average_rad,2)*E;
					electricproperties->chargeM *= (1.0 - 0.3*ChRandom() );
					
					}
				}

				ChVector<> ElectricForce = electricproperties->ElectricForce;

				electricproperties->ElectricForce = 0;

				// switch off electric forces if too out-of-plane
				if ((mrelpos.z > conveyor_width*0.5) || (mrelpos.z < -conveyor_width*0.5))
					ElectricForce = 0; 

				abody->Accumulate_force(ElectricForce, abody->GetPos(), false);


			} // end if material==metal


			//
			//===== PLASTIC FORCES ==========
			//

		    

			if (electricproperties->material_type == ElectricParticleProperty::e_mat_plastic) //forze sulle particelle non metalliche
			{
				

				// charge the particle? (contact w. drum)
				if ((distx > 0.04) && (disty > 0))
				{
					if (electricproperties->chargeP == 0)
					{
						electricproperties->chargeP = 3*CH_C_PI*epsilonO*pow(2*average_rad,2)*450000*(epsilonR/(epsilonR+2)); // charge
						electricproperties->chargeP *= (1.0 - 0.3*ChRandom() );
					}
				} 
				// discharge the particle? (contact w. blade)
				if (distx < -(drumdiameter*0.5 -0.009) && (disty > -(drumdiameter*0.5 + 0.009)) || sqrt(pow(distx,2)+ pow(disty,2))> (1.03*drumdiameter*0.5))
				{
					electricproperties->chargeP = 0; // charge
				}

				ChVector<> ElectricImageForce = electricproperties->ElectricImageForce;


				electricproperties->ElectricImageForce.x = 0;
				electricproperties->ElectricImageForce.y = 0;
				electricproperties->ElectricImageForce.z = 0;	
						
				
				// switch off electric forces if too out-of-plane
				if ((mrelpos.z > conveyor_width*0.5) || (mrelpos.z < -conveyor_width*0.5))
					ElectricImageForce = 0; 


				abody->Accumulate_force(ElectricImageForce, abody->GetPos(), false);
				

			   }  // end if material==plastic
			




			//ChVector<> DragForce;
			//DragForce.x = -CD*ro*velocity_norm_sq*CH_C_PI*diam.x*diam.y/2*cos(phi2);
			//DragForce.y = CD*ro*velocity_norm_sq*CH_C_PI*diam.x*diam.y/2*sin(phi2);
			//DragForce.z = 0;

			//abody->Accumulate_force( DragForce, abody->GetPos(), false);

			//ChVector<> LiftForce;
			//LiftForce.x = CL*ro*velocity_norm_sq*CH_C_PI*diam.x*diam.y/2*sin(phi2);
			//LiftForce.y = CL*ro*velocity_norm_sq*CH_C_PI*diam.x*diam.y/2*cos(phi2);
			//LiftForce.z = 0;	
		
			//abody->Accumulate_force(LiftForce, abody->GetPos(), false);
			

		} // end if(was_a_particle) , i.e. a body with electrical asset

	} // end for() loop on all bodies
}
*/ 


void draw_forces(ChIrrApp& application, double scalefactor = 1.0)
{
	ChSystem* msystem = application.GetSystem();

	for (unsigned int i=0; i<msystem->Get_bodylist()->size(); i++)
	{
		ChBody* abody = (*msystem->Get_bodylist())[i];

		bool was_a_particle = false;
		ChSharedPtr<ElectricParticleProperty> electricproperties; // null by default

		// Fetch the ElectricParticleProperty asset from the list of 
		// assets that have been attached to the object, and retrieve the
		// custom data that have been stored. ***ALEX
		for (unsigned int na= 0; na< abody->GetAssets().size(); na++)
		{
			ChSharedPtr<ChAsset> myasset = abody->GetAssetN(na);
			if (myasset.IsType<ElectricParticleProperty>())
			{
				// OK! THIS WAS A PARTICLE! ***ALEX
				was_a_particle = true;		
				electricproperties = myasset.DynamicCastTo<ElectricParticleProperty>();
			} 
		}

		// Do the computation of forces only on bodies that had 
		// the 'ElectricParticleProperty' attached.. **ALEX
		if(was_a_particle)
		{
			ChVector<> custom_force = abody->Get_accumulated_force();
			custom_force *= scalefactor;
			ChIrrTools::drawSegment(application.GetVideoDriver(), 
				abody->GetPos(), 
				abody->GetPos() + custom_force,
				video::SColor(255,   0,0,255));
		}
	}
}


void UpdateTrajectories(ChIrrApp& application)
{
	ChSystem* msystem = application.GetSystem();

	for (unsigned int i=0; i<msystem->Get_bodylist()->size(); i++)
	{
		ChBody* abody = (*msystem->Get_bodylist())[i];

		bool was_a_particle = false;
		ChSharedPtr<ParticleTrajectory> trajectoryasset; // null by default

		// Fetch the ElectricParticleProperty asset from the list of 
		// assets that have been attached to the object, and retrieve the
		// custom data that have been stored. ***ALEX
		for (unsigned int na= 0; na< abody->GetAssets().size(); na++)
		{
			ChSharedPtr<ChAsset> myasset = abody->GetAssetN(na);
			if (myasset.IsType<ParticleTrajectory>())
			{
				// OK! trajectory storage!	
				trajectoryasset = myasset.DynamicCastTo<ParticleTrajectory>();
				trajectoryasset->positions.push_back( abody->GetPos() );
				trajectoryasset->speeds.push_back( abody->GetPos_dt() );
				
				// remove excessive amount of elements
				while (trajectoryasset->positions.size() > trajectoryasset->max_points )
					trajectoryasset->positions.pop_front();
				while (trajectoryasset->speeds.size() > trajectoryasset->max_points )
					trajectoryasset->speeds.pop_front();
			} 
		}
	}
}

void DrawTrajectories(ChIrrApp& application)
{
	ChSystem* msystem = application.GetSystem();

	for (unsigned int i=0; i<msystem->Get_bodylist()->size(); i++)
	{
		ChBody* abody = (*msystem->Get_bodylist())[i];

		ChVector<> pointA = abody->GetPos();

		bool was_a_particle = false;
		ChSharedPtr<ParticleTrajectory> trajectoryasset; // null by default

		// Fetch the ElectricParticleProperty asset from the list of 
		// assets that have been attached to the object, and retrieve the
		// custom data that have been stored. ***ALEX
		for (unsigned int na= 0; na< abody->GetAssets().size(); na++)
		{
			ChSharedPtr<ChAsset> myasset = abody->GetAssetN(na);
			if (myasset.IsType<ParticleTrajectory>())
			{
				trajectoryasset = myasset.DynamicCastTo<ParticleTrajectory>();
				int npoints = 0;
				std::list< ChVector<> >::const_iterator iterator;
				std::list< ChVector<> >::const_iterator iteratorspeed;
				iteratorspeed = trajectoryasset->speeds.begin();
				for (iterator = trajectoryasset->positions.begin(); iterator != trajectoryasset->positions.end(); ++iterator)
				{
					ChVector<> pointB = *iterator;
					ChVector<> speed = *iteratorspeed;
					if (npoints >0)
					{
						double scalarspeed = speed.Length();
						double normalizedspeed = scalarspeed / 5.0;
						video::SColor mcol (255, (int)(255.*normalizedspeed), (int)(255.*normalizedspeed), (int)(255.*(1.0-normalizedspeed)) );
						ChIrrTools::drawSegment(application.GetVideoDriver(), 
								pointA, 
								pointB,
								mcol);
					}
					pointA = pointB;
					++npoints;
					++iteratorspeed;
				}
			} 
		}
	}
}



int main(int argc, char* argv[])
{
	// In CHRONO engine, The DLL_CreateGlobals() - DLL_DeleteGlobals(); pair is needed if
	// global functions are needed. 
	DLL_CreateGlobals();

	// Create a ChronoENGINE physical system
	ChSystem mphysicalSystem;

	// From now on, functions in ChParticlesSceneNodeTools will find 3d .obj models 
	// in "../objects/", instead of default "../data/" dir:
	irrlicht_default_obj_dir = "../objects/";

	// Create the Irrlicht visualization (open the Irrlicht device, 
	// bind a simple user interface, etc. etc.)
	ChIrrApp application(&mphysicalSystem, L"Conveyor belt",core::dimension2d<u32>(800,600),false);

	// Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
	//application.AddTypicalLogo("../objects/");
	application.AddTypicalSky("../objects/skybox/");
	application.AddTypicalLights();
	application.AddTypicalCamera(core::vector3df(1.5f,0.4f,-1.0f), core::vector3df(0.5f,0.f,0.f));
	if (irr_cast_shadows)
		application.AddLightWithShadow(vector3df(-4.5f,5.5f,4.5f), vector3df(0.f,0.f,0.f), 10, 1.2,10.2, 30,512, video::SColorf(1.f,0.9f,0.9f));

	// This is for GUI tweaking of system parameters..
	MyEventReceiver receiver(&application);
	  // note how to add the custom event receiver to the default interface:
	application.SetUserEventReceiver(&receiver);


	// Set small collision envelopes for objects that will be created from now on..
	ChCollisionModel::SetDefaultSuggestedEnvelope(0.002); 
	//ChCollisionModel::SetDefaultSuggestedMargin  (0.002);
	ChCollisionModel::SetDefaultSuggestedMargin  (0.0008);



	

	// IMPORT A SOLIDWORK MODEL 

	// 1) create the Python engine. This is necessary in order to parse the files that 
	// have been saved using the SolidWorks add-in for Chrono::Engine.

	ChPythonEngine my_python;

	// 2) loads the .py file (as saved from SolidWorks) and fill the system.
	try
	{
		//my_python.ImportSolidWorksSystem("../CAD_conveyor/conveyor_Ida", mphysicalSystem);  // note, don't type the .py suffix in filename..
		//my_python.ImportSolidWorksSystem("../CAD_conveyor/prova_ida_senzacinghia", mphysicalSystem);
		//my_python.ImportSolidWorksSystem("../CAD_conveyor/prova_edo", mphysicalSystem);
		//my_python.ImportSolidWorksSystem("../CAD_conveyor/edo_mod_1206", mphysicalSystem);
		my_python.ImportSolidWorksSystem("../CAD_conveyor/Edo_0707_limits", mphysicalSystem);
	}
	catch (ChException myerror)
	{
		GetLog() << myerror.what();
	}


	// 3) fetch coordinate values and objects from what was imported from CAD

	ChCoordsys<> conveyor_csys = CSYSNORM;


	ChSharedPtr<ChMarker> my_marker = mphysicalSystem.SearchMarker("conveyor_origin").DynamicCastTo<ChMarker>();
	if (my_marker.IsNull())
		GetLog() << "Error: cannot find conveyor_origin marker from its name in the C::E system! \n";
	else
		conveyor_csys = my_marker->GetAbsCoord();  // fetch both pos and rotation of CAD
		
	
	my_marker = mphysicalSystem.SearchMarker("Particles_origin_new").DynamicCastTo<ChMarker>();
	if (my_marker.IsNull())
		GetLog() << "Error: cannot find particles_origin marker from its name in the C::E system! \n";
	else
		nozzle_csys = my_marker->GetAbsCoord();  // fetch both pos and rotation of CAD

    my_marker = mphysicalSystem.SearchMarker("Center_rotor").DynamicCastTo<ChMarker>();
	if (my_marker.IsNull())
		GetLog() << "Error: cannot find centro_cilindro marker from its name in the C::E system! \n";
	else
		drum_csys = my_marker->GetAbsCoord();  // fetch both pos and rotation of CAD

		// fetch mrigidBodyDrum pointer! will be used for changing the friction, the collision family, and later to create the motor
	ChSharedPtr<ChBodyAuxRef> mrigidBodyDrum = mphysicalSystem.Search("Slave-1/hub_big_slave-1").DynamicCastTo<ChBodyAuxRef>();  
	if (mrigidBodyDrum.IsNull())
		GetLog() << "ERROR: cannot find drum from its name in the C::E system! ! \n";
	else
	{
		mrigidBodyDrum->GetCollisionModel()->SetFamily(7);
		mrigidBodyDrum->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(8);
		mrigidBodyDrum->SetFriction(0.2f); 
		//mrigidBodyDrum->SetImpactC(0.0f);
		//mrigidBodyDrum->SetRollingFriction(0.2f);
		//mrigidBodyDrum->SetSpinningFriction(0.2f);
	}
	

	//
	// Create a truss (absolute fixed reference body, for connecting the rotating cyl.)
	//

	ChSharedPtr<ChBody> mtruss(new ChBody);
	mtruss->SetBodyFixed(true);

	// Finally, do not forget to add the body to the system:
	application.GetSystem()->Add(mtruss);
   

	//
	// Create the conveyor belt (this is a pure Chrono::Engine object, 
	// because an Irrlicht 'SceneNode' wrapper is not yet available, so it is invisible - no 3D preview)
	
	// NOTE!!! for the Corona Effect Separator this should NOT be needed....

	ChSharedPtr<ChConveyor> mconveyor;

	mconveyor = ChSharedPtr<ChConveyor>(new ChConveyor(conveyor_length, conv_thick, conveyor_width));
	mconveyor->SetBodyFixed(true);
	mconveyor->SetFriction(0.6f);
	//mconveyor->SetRollingFriction(0.00f);
	mconveyor->SetConveyorSpeed(STATIC_speed);
	mconveyor->SetCollide(true);

	mconveyor->SetPos( conveyor_csys.pos ); // mconveyor->SetCoord ( conveyor_csys );
	mconveyor->ConcatenatePreTransformation(ChFrameMoving<>(ChVector<>(0, -conv_thick/2., 0))); // cut if csys already in the middle 

	//mconveyor->GetPlate()->GetCollisionModel()->SetFamily(8); // note the additional  ->GetPlate()
	mconveyor->GetPlate()->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(7); 

	mphysicalSystem.Add(mconveyor);

	// Attach a visualization shape asset (optional). // Not needed - already in .obj meshes from SW
	
	ChSharedPtr<ChBoxShape> mbox(new ChBoxShape);
	mbox->GetBoxGeometry().SetLenghts(ChVector<>(conveyor_length,conv_thick, conveyor_width));
	mconveyor->AddAsset(mbox);
	

	// //***TEST*** mistero!!! (spiegato) si può togliere...
	//  ChBodySceneNode* mfence2 = (ChBodySceneNode*)addChBodySceneNode_easyBox(
	//											&mphysicalSystem, application.GetSceneManager(),
	//											1.0,
	//											conveyor_csys.pos + ChVector<>(0, 0.2, 0),
	//											ChQuaternion<>(1,0,0,0), 
	//											ChVector<>(0.1,0.11,0.1) );
	// //mfence2->GetBody()->SetBodyFixed(true);
	// mfence2->GetBody()->SetFriction(0.1);


	// 
	// Create a motor constraint between the cylinder and the truss
	//

/*  ***TODO*** reactivate this after drum_csys has a reasonable value from SW.... 
	ChSharedPtr<ChLinkEngine> mengine;

	if (!mrigidBodyDrum.IsNull())
	{
		mengine = ChSharedPtr<ChLinkEngine>(new ChLinkEngine);
		ChSharedPtr<ChBody> mdrum(mrigidBodyDrum);
		mengine->Initialize(mdrum, mtruss, drum_csys);

		mengine->Set_eng_mode(ChLinkEngine::ENG_MODE_SPEED);
		if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(mengine->Get_spe_funct()))
			mfun->Set_yconst(-STATIC_speed/(drumdiameter*0.5));

		// Finally, do not forget to add the body to the system:
		application.GetSystem()->Add(mengine);
	}
*/

  



	//
	// Create an (optional) exporter to POVray 
	// 

	ChPovRay pov_exporter = ChPovRay(&mphysicalSystem);

	if (save_POV_screenshots)
	{
				// Sets some file names for in-out processes.
		pov_exporter.SetTemplateFile		("../objects/_template_POV.pov");
		pov_exporter.SetOutputScriptFile    ("rendering_frames.pov");
			
				// save the .dat files and the .bmp files
				// in two subdirectories, to avoid cluttering the current directory...
		ChFileutils::MakeDirectory("outputPOV");
		ChFileutils::MakeDirectory("animPOV");

		pov_exporter.SetOutputDataFilebase("outputPOV/my_state");
		pov_exporter.SetPictureFilebase("animPOV/picture");

				// optional: modify the POV default light
		pov_exporter.SetLight(ChVector<>(1.5f,4.4f,-1.0f), ChColor(0.1f,0.1f,0.1f), false);

		pov_exporter.SetCamera(ChVector<>(1.5f,0.8f,-1.0f),ChVector<>(0.5f,0.f,0.f),60,false);

				// optional: use SetCustomPOVcommandsScript() to add further POV commands,
				// ex. create an additional light, and an additional grid, etc. 
				// Remember the "\" char per each newline.
		
		pov_exporter.SetCustomPOVcommandsScript(" \
			light_source {   \
			  <2, 10, -3>  \
			  color rgb<1.8,1.8,1.8> \
			  area_light <4, 0, 0>, <0, 0, 4>, 5, 5 \
			  adaptive 1 \
			  jitter\
			} \
			object{ Grid(0.5,0.01, rgb<0.9,0.9,0.9>, rgbt<1,1,1,1>) rotate <90, 0, 0>  } \
			object{ Grid(0.1,0.04, rgb<1.5,1.5,1.5>, rgbt<1,1,1,1>) rotate <90, 0, 0> translate 0.001*z} \
		");

				// IMPORTANT! Tell to the POVray exporter that 
				// he must take care of converting the shapes of
				// all items!
		pov_exporter.AddAll();

				// IMPORTANT! Create the two .pov and .ini files for POV-Ray (this must be done
				// only once at the beginning of the simulation).
		pov_exporter.ExportScript();

	}


	//
	// For enabling Irrlicht visualization of assets (that have been added so far)
	//

	application.AssetBindAll();
	application.AssetUpdateAll();
	if (irr_cast_shadows)
		application.AddShadowAll();


	// 
	// THE SOFT-REAL-TIME CYCLE
	//
	
	application.SetStepManage(true);
	application.SetTimestep(0.001);
	
	application.GetSystem()->SetIntegrationType(ChSystem::INT_ANITESCU);
	application.GetSystem()->SetLcpSolverType(ChSystem::LCP_ITERATIVE_SOR_MULTITHREAD); // or ChSystem::LCP_ITERATIVE_BARZILAIBORWEIN for max precision
		// important! dt is small, and particles are small, so it's better to keep this small...
	application.GetSystem()->SetMaxPenetrationRecoverySpeed(0.15);// not needed in INT_TASORA, only for INT_ANITESCU
	application.GetSystem()->SetMinBounceSpeed(0.01);

	application.GetSystem()->Set_G_acc(ChVector<>(0, -9.81, 0));

	double threshold = -0.815;
	ofstream record;
	record.open("output\\threshold.dat",ios::trunc); 
	record << "Time" << "\t";
	record << "X Position" << "\t";
	record << "Y Position" << "\t";
	record << "Z Position" << "\t";
	record << "Shape" << "\t";
	record << "Mass (kg)" << "\t";
	record << "Conductivity" << "\n";
	
	int savenum = 0;

	ChFileutils::MakeDirectory("screenshots");
	ChFileutils::MakeDirectory("output");
	
	application.GetSystem()->ShowHierarchy(GetLog());

	while(application.GetDevice()->run()) 
	{
		application.GetVideoDriver()->beginScene(true, true, SColor(255,140,161,192));

		application.DrawAll();

		application.DoStep();

		if (!application.GetPaused())
		{ 
			
			totframes++;

			apply_forces (	&mphysicalSystem,		// contains all bodies
							drum_csys,		 // pos and rotation of axis of drum (not rotating reference!)
							drumspeed,		 // speed of drum
						    numberofpoles,   // ***Edo
							intensity,       // ***Edo
							drumdiameter,    // ***Edo
							totframes);      // ***Edo
							//h1,            // useless for ECS (at least for the moment) ***Edo
							//h2,
							//L,
							//electrodediameter,
							//j,
							//alpha,
							//U,
							//f,
							

		

			if (receiver.checkbox_plotforces->isChecked())
				draw_forces ( application , 1000);

			if (receiver.checkbox_plottrajectories->isChecked())
				DrawTrajectories(application);

			
			// Continuosly create debris that fall on the conveyor belt
			
			if (debris_number <= max_numb_particles)
			{
			create_debris(	application.GetTimestep(), 
							STATIC_flow, 
							nozzle_csys,
							*application.GetSystem(), 
							&application, 
							&pov_exporter);

			debris_number = debris_number + particles_dt;
			}

			// Limit the max age (in seconds) of debris particles on the scene, 
			// deleting the oldest ones, for performance
			purge_debris (*application.GetSystem(),6);

			// Maybe the user played with the slider and changed STATIC_speed...
			if (!mconveyor.IsNull())
		 		mconveyor->SetConveyorSpeed(STATIC_speed);

			// ..also for rotating drum:
			/*  ***TODO*** reactivate this after drum_csys has a reasonable value from SW.... 
			if (!mengine.IsNull())
			  if (ChFunction_Const* mfun = dynamic_cast<ChFunction_Const*>(mengine->Get_spe_funct()))
				mfun->Set_yconst(-STATIC_speed/(drumdiameter*0.5));
			*/

			// update the assets containing the trajectories, if any
			if (receiver.checkbox_plottrajectories->isChecked())
				if (totframes % 20 == 0)
					UpdateTrajectories(application);

			// Save data on file (each n integration steps, to avoid filling
			// the hard disk and to improve performance)
			if (totframes % saveEachNframes == 0)
			{
				savenum++;

				// Save log file as '.txt' files?

				if (save_dataset == true)
				{
					char buffer[120];
					sprintf(buffer, "output/esempio_output%05d.txt", savenum);
					GetLog() << "\n saving dataset: " << buffer;
					ChStreamOutAsciiFile file_for_output(buffer);
					for (unsigned int i=0; i<mphysicalSystem.Get_bodylist()->size(); i++)
					{
						ChBody* abody = (*mphysicalSystem.Get_bodylist())[i];

						// Fetch the ElectricParticleProperty asset from the list
						for (unsigned int na= 0; na< abody->GetAssets().size(); na++)
						{
							ChSharedPtr<ChAsset> myasset = abody->GetAssetN(na);
							
							if (myasset.IsType<ElectricParticleProperty>())
							{
								// ok, its a particle!

								
								ChSharedPtr<ElectricParticleProperty> electricproperties = myasset.DynamicCastTo<ElectricParticleProperty>();
								//double my_cond  = electricproperties->conductivity ;
								ChVector<> my_InducedForce = electricproperties->InducedForce;
								ChVector<> my_InducedTorque = electricproperties->InducedTorque;
								ChVector<> my_DragForce = electricproperties->DragForce;
								ChVector<> my_LiftForce = electricproperties->LiftForce;
								
								double rad = ((abody->GetMass())*3)/((abody->GetDensity())*4*CH_C_PI);
								
								// Save on disk some infos...
								file_for_output << abody->GetIdentifier()<< ", "
												<< abody->GetPos().x << ", "
												<< abody->GetPos().y << ", "
												<< abody->GetPos().z << ", "
												<< abody->GetDensity() <<"\n";
												//<< my_cond << ", "
												//<< abody->GetMass()<< ", "
												//<< pow(rad,1.0/3) << "\n";
                                                //<< abody->GetPos_dt().x << ", "
												//<< abody->GetPos_dt().y << ", "
												//<< abody->GetPos_dt().z << ", "
											    //<< my_InducedForce << ", "
												//<< my_InducedTorque << ", "
                                                //<< my_DragForce <<  ", " 
												//<< my_LiftForce << "\n";
						                       
								               
							}
						}
					}
				}
				
				// Save Irrlicht screenshots?

				if (save_irrlicht_screenshots == true)
				{
					video::IImage* image = application.GetVideoDriver()->createScreenShot();
					char buffer[120];
					sprintf(buffer, "screenshots/screenshot%05d.bmp", savenum);
					GetLog() << "\n saving screenshot: " << buffer;
					if (image)
						application.GetVideoDriver()->writeImageToFile(image, buffer);
					image->drop();
				}

				// Save POV screenshots?

				if(save_POV_screenshots)
				{
					pov_exporter.ExportData();
					GetLog() << "\n saving POV data n." << savenum;
				}

			} // end saving code

			
		
		}

		application.GetVideoDriver()->endScene();  
		
	}
	

	

 
	// Remember this at the end of the program, if you started
	// with DLL_CreateGlobals();
	DLL_DeleteGlobals();

	return 0;
}