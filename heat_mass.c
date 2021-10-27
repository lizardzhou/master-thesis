/***********************************************************************
   UDF for defining the heat and mass transport for
   multicomponent particle vaporization
 ***********************************************************************/
#include "udf.h"
#include "dpm_mem.h"
#define x_cr 0.7 /* critical moisture content, end of 1st drying stage */
#define x_eq 0.05 /* final moisture content, end of 2nd drying stage */
#define kappa 2 /* curvature of normalized drying curve */

DEFINE_INIT(dpm_setup, domain)
{
	/* Allocate the memory for self-defined variables */
	if (NULLP(user_particle_vars)) {
		Init_User_Particle_Vars();
	}
	/* Set the name of self-defined scalar variables: moisture content, relative drying rate, particle diameter */
	strcpy(user_particle_vars[0].name, "moisture-content");
	strcpy(user_particle_vars[0].label, "Moisture Content");
	strcpy(user_particle_vars[1].name, "relative-drying-rate");
	strcpy(user_particle_vars[1].label, "Relative Drying Rate");
	strcpy(user_particle_vars[2].name, "particle-diameter");
	strcpy(user_particle_vars[2].label, "Particle Diameter from CDC"); /* CDC: characteristic drying curve */
}

DEFINE_DPM_SCALAR_UPDATE(new_variables, c, t, init, tp)
{
	int ns;
	int nc = TP_N_COMPONENTS(tp);	/*number of components in particle*/
	real xs = 0.21; /* initial solid mass fraction in particle */
	real xw0 = 1 - xs; /* initial water mass fraction in particle */
	real mp0 = TP_INIT_MASS(tp); /* initial particle mass */
	real mp = TP_MASS(tp); /* particle mass */
	real xv, wap_mass, xi, f, dp; /* xv: water mass fraction in particle. wap_mass: water mass in particle. xi: moisture content in particle. f: relative drying rate. dp: particle diameter */
	
	if (init) /* initial values */
	{
		/* moisture content in particle */
		TP_USER_REAL(tp, 0) = mp0 * xw0 / (mp0 - mp0 * xw0);
		/* relative drying rate */
		TP_USER_REAL(tp, 1) = 1;
		/* particle diameter according to CDC */
		TP_USER_REAL(tp, 2) = TP_INIT_DIAM(tp);
	}
	else /* values during simulation */
	{
		for (ns = 0; ns < nc; ns++) /* loop for every vaporating species */
		{
			int gas_index = TP_COMPONENT_INDEX_I(tp, ns); /*find the index of evaporating species*/
			if (ns == gas_index) /*IMPORTANT: this part works ONLY WHEN the water component is set to 1st in injection "multicomponent" option!!! */
			{
				xv = TP_COMPONENT_I(tp, ns); /*mass fraction of evaporated species in particle*/
				wap_mass = mp * xv; /*water mass in particle*/
				xi = wap_mass / (mp - wap_mass); /*moisture content in particle*/

				if (xi >= x_cr) /*1st drying stage*/
				{
					f = 1.0;
					dp = TP_DIAM(tp);
				}
				else if (xi > x_eq) /*2nd drying stage*/
				{
					real mu = (xi - x_eq) / (x_cr - x_eq);
					f = kappa * mu / (1 + (kappa - 1) * mu) ;
					dp = TP_USER_REAL(tp, 2);
				}
				else /*3rd drying stage*/
				{
					xi = x_eq;
					f = 0;
					dp = TP_USER_REAL(tp, 2);
				}
			}
		}

		/* assign values after calculation */
		TP_USER_REAL(tp, 0) = xi; /* calculated moisture content*/
		TP_USER_REAL(tp, 1) = f; /* calculated relative drying rate */
		TP_USER_REAL(tp, 2) = dp; /* calculated particle diameter */
	}
}


DEFINE_DPM_HEAT_MASS(multivap,tp,Cp,hgas,hvap,cvap_surf,Z,dydt,dzdt)
{
	int ns;
	Material *sp;
	real dens_total = 0.0;     /* total vapor density*/
	real P_total = 0.0;      /* vapor pressure */
	int nc = TP_N_COMPONENTS(tp);   /* number of particle components */
	Thread *t0 = TP_CELL_THREAD(tp);   /* thread where the particle is in*/
	
	Material *gas_mix = THREAD_MATERIAL(DPM_THREAD(t0, tp)); /* gas mixture material */
	Material *cond_mix = TP_MATERIAL(tp); /* particle mixture material*/
	cphase_state_t *c = &(tp->cphase[0]); /* cell information of particle location*/
	real molwt[MAX_SPE_EQNS]; /* molecular weight of gas species */
	real Tp = TP_T(tp);   /* particle temperature */
	real mp = TP_MASS(tp);   /* particle mass */
	real molwt_bulk = 0.;  /* average molecular weight in bulk gas */
	real Dp = DPM_DIAM_FROM_VOL(mp / TP_RHO(tp)); /* particle diameter */
	real dpp = TP_USER_REAL(tp, 2); /*particle diameter (user-defined)*/
	real Ap = DPM_AREA(dpp);      /* particle surface area (user-defined)*/
	real f = TP_USER_REAL(tp, 1); /* relative drying rate (user-defined)*/
	real Pr = c->sHeat * c->mu / c->tCond;   /* Prandtl number (heat transfer) = Cp*mu/lambda */
	real Nu = 2.0 + 0.664 * sqrt(tp->Re * (dpp / Dp)) * pow(Pr, 1. / 3.); /* Nusselt number (heat transfer) */
	real h = Nu * c->tCond / dpp;     /* Heat transfer coefficient = Nu*lambda/characteristicLength(particleDiameter) [W/(m2*K)])*/
	real dh_dt = h * (c->temp - Tp) * Ap;  /* heat source term [W]*/
	real xv, wap_mass, xi; /* xv: water mass fraction in particle. wap_mass: water mass in particle. xi: moisture content in particle. (user-defined) */
	dydt[0] += dh_dt / (mp * Cp); /* particle temperature [K/s] */
	dzdt->energy -= dh_dt; /* gas phase enthalpy [W] */
	
	mixture_species_loop(gas_mix,sp,ns)
	{
		molwt[ns] = MATERIAL_PROP(sp,PROP_mwi);/* molecular weight of gas species */
		molwt_bulk += c->yi[ns] / molwt[ns];/* average molecular weight */
	}
	
	/* prevent division by zero */
	molwt_bulk = MAX(molwt_bulk,DPM_SMALL);

	
	for (ns = 0; ns < nc; ns++) /* loop for every vaporating species */
	{
		int gas_index = TP_COMPONENT_INDEX_I(tp,ns);  /* gas species index of vaporization */
		/* water mass fraction in particle */
		xv = TP_COMPONENT_I(tp,ns);
		/* water mass in particle */
		wap_mass = mp * xv;
		/* moisture content in particle */
		xi = wap_mass / (mp - wap_mass); 
		
		if(gas_index >= 0)
		{
			/* condensed material */
			Material *cond_c = MIXTURE_COMPONENT(cond_mix, ns);
			/* vaporization temperature */
			real vap_temp = MATERIAL_PROP(cond_c,PROP_vap_temp);
			/* diffusion coefficient [m2/s] */
			real D = DPM_BINARY_DIFFUSIVITY(tp,cond_c,TP_T(tp));
			/* Schmidt number (mass transfer) = mu/(rho*difusionCoefficient) */
			real Sc = c->mu / (c->rho * D);
			/* Sherwood number (mass transfer) */
			real Sh = 2. + 0.664 * sqrt((tp->Re) * (dpp / Dp)) * pow(Sc, 1. / 3.);
			/* mass transfer coefficient = Sh*diffusionCoefficient/characteristicLength(particleDiameter) [m/s]*/
			real k = Sh * D / dpp;
			/* bulk gas concentration (ideal gas) [mol/m3] */
			real cvap_bulk = c->pressure / UNIVERSAL_GAS_CONSTANT / c->temp * c->yi[gas_index] / molwt_bulk / solver_par.molWeight[gas_index];
			/* vaporization rate [kg/s] */
			real vap_rate = f * k * molwt[gas_index] * Ap * (cvap_surf[ns] - cvap_bulk);
			/* no vaporization below vaporization temperature, and no condensation */
			if (Tp < vap_temp || vap_rate < 0.0)
			{
				vap_rate = 0.;
			}

			/* particle component mass [kg/s] */
			dydt[1+ns] -= vap_rate;
			/* gas phase species mass [kg/s] */
			dzdt->species[gas_index] += vap_rate;
			/* dT/dt = dh/dt / (m*Cp) [K/s] */
			dydt[0] -= hvap[gas_index] * vap_rate / (mp * Cp);
			/* gas enthalpy source term [W] */
			dzdt->energy += hgas[gas_index] * vap_rate;

			P_total += cvap_surf[ns]; /* update vapor pressure */
			dens_total += cvap_surf[ns] * molwt[gas_index]; /* update total vapor density */
		}
	}
	
}

