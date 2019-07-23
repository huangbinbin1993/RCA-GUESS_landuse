///////////////////////////////////////////////////////////////////////////////////////
/// \file landcover.cpp
/// \brief Functions handling landcover aspects, such as creating or resizing Stands
///
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "landcover.h"
#include "guessio.h"

void landcover_init(Gridcell& gridcell) {
	landcovertype landcover;

	getlandcover(gridcell);		//Gets gridcell.landcoverfrac from landcover input file(s) or ins-file.

	for(int i=0;i<NLANDCOVERTYPES;i++) { //For all landcover types without subclasses
//		if(i!=CROPLAND) {					// cropland subclasses turned off in this version
			if(run[i]) {
				// We always create the stands in RCA, regardless of whether they have
				// a non-zero fraction. There is code that assumes we always have a forest
				// and an open land stand (for instance when reading in climate).
//				if(gridcell.landcoverfrac[i]>0.0) {
					landcover=(landcovertype)i;
					Stand& stand=gridcell.createobj(gridcell,landcover);

					pftlist.firstobj();
					while (pftlist.isobj) {
						Pft& pft=pftlist.getobj();
						if(pft.landcover==i) {
							stand.pft[pft.id].active=true;
						}
						pftlist.nextobj();
					}
//				}
			}
//		}
	}
}

void harvest_natural(double& cmass_leaf,double& cmass_root,double& cmass_sap,double& cmass_heart,double& cmass_debt,
	double& litter_leaf,double& litter_root,double& litter_sap,double& litter_heart,double& acflux_harvest,double& harvested_products_slow,Individual& indiv) 
{
	double harvest=0.0;
	double residue_outtake=0.0;
	bool alive=indiv.alive;

	if(alive && cmass_root>0.0)
		litter_root+=cmass_root;			//all root carbon goes to litter
	cmass_root=0.0;

	if(alive && (cmass_sap+cmass_heart-cmass_debt)>0.0)						// Only wood currently harvested in this function !
	{	//indiv.pft.harv_eff = Bondeau's removal = 0.7 for tree wood, 0 for the rest
		harvest=indiv.pft.harv_eff*(cmass_sap+cmass_heart-cmass_debt);		//harvested products

		if(ifslowharvestpool)
		{
			harvested_products_slow+=harvest*indiv.pft.harvest_slow_frac;	//harvested products not consumed (oxidized) this year put into patchpft.harvested_products_slow
			harvest=harvest*(1-indiv.pft.harvest_slow_frac);
		}

		acflux_harvest+=harvest;							//harvested products consumed (oxidized) this year put into patch.fluxes.acflux_harvest, not litter pool !

		cmass_sap=(1-indiv.pft.harv_eff)*cmass_sap;			//unharvested parts of the plant
		cmass_heart=(1-indiv.pft.harv_eff)*cmass_heart;
		cmass_debt=(1-indiv.pft.harv_eff)*cmass_debt;		// ????
	}

	residue_outtake=indiv.pft.res_outtake*(cmass_sap+cmass_heart-cmass_debt+cmass_leaf);
	acflux_harvest+=residue_outtake;																//removed residues

	litter_leaf+=cmass_leaf*(1-indiv.pft.res_outtake);												//not removed residues
	litter_sap+=cmass_sap*(1-indiv.pft.res_outtake);												//not removed residues
	litter_heart+=(cmass_heart-cmass_debt)*(1-indiv.pft.res_outtake);								//not removed residues

	cmass_sap=cmass_heart=cmass_debt=cmass_leaf=0.0;
}

void landcover_dynamics(Gridcell& gridcell)
{	// Called first day of the year if run_landcover is set.
	int i;	
	landcovertype landcover;
	double landcoverfrac_change[NLANDCOVERTYPES]={0.0};
//	double cropfrac_change[NCROPSTANDS_MAX]={0.0};	
	double cropfrac_sum_old=0.0;
//	double cropstand_change[NCROPSTANDS_MAX]={0.0};

	gridcell.LC_updated=false;

/////////////////////////////////////////////////////////////
//Landcover and cft fraction update (from updated landcoverfrac):

//Save old fraction values:
	for(i=0;i<NLANDCOVERTYPES;i++)
		gridcell.landcoverfrac_old[i]=gridcell.landcoverfrac[i];
//	for(i=0;i<NCROPSTANDS_MAX;i++)
//		cropfrac_sum_old+=gridcell.cftfrac_old[i]=gridcell.cftfrac[i];

//Get new gridcell.landcoverfrac and/or gridcell.cftfrac from LUdata and CFTdata.
	if(!all_fracs_const)
		getlandcover(gridcell);	
	else return;

	// Don't do the rest of landcover_dynamics in RCA.
	// This version of GUESS (pre-crop) doesn't do it properly yet anyway,
	// and we don't want to remove stands if their fraction becomes zero
	// (because we have assumed in other places that there are always two
	// stands, this assumption can perhaps be removed in future)
	return;

	double changeLC=0.0;
	double change_crop=0.0;
	double change_stand=0.0;
	double transferred_fraction=0.0;
	double receiving_fraction=0.0;

	if(!lcfrac_fixed)
	{
		for(i=0;i<NLANDCOVERTYPES;i++)
		{
			landcoverfrac_change[i]=gridcell.landcoverfrac[i]-gridcell.landcoverfrac_old[i];
			changeLC+=fabs(landcoverfrac_change[i])/2.0;
//			if(i!=CROPLAND)											//Landcovers with only one stand.
			{
				if(landcoverfrac_change[i]<0.0)
					transferred_fraction-=landcoverfrac_change[i];
				if(landcoverfrac_change[i]>0.0)
					receiving_fraction+=landcoverfrac_change[i];
				change_stand+=fabs(landcoverfrac_change[i])/2.0;
			}
		}
	}
/*	if(run[CROPLAND] &&!cftfrac_fixed)
	{
		for(i=0;i<NCROPSTANDS_MAX;i++)
		{
			cropfrac_change[i]=gridcell.cftfrac[i]-gridcell.cftfrac_old[i];
			cropstand_change[i]=gridcell.cftfrac[i]*gridcell.landcoverfrac[CROPLAND]-gridcell.cftfrac_old[i]*gridcell.landcoverfrac_old[CROPLAND];

			if(cropstand_change[i]<0.0)
				transferred_fraction-=cropstand_change[i];
			if(cropstand_change[i]>0.0)
				receiving_fraction+=cropstand_change[i];

			if(cropfrac_sum_old!=0.0)
			{
				change_crop+=fabs(cropfrac_change[i])/2.0;
			}
			else
			{
				change_crop+=fabs(cropfrac_change[i]);
			}
			change_stand+=fabs(cropstand_change[i])/2.0;	//cropfrac_sum_old+gridcell.landcoverfrac[NATURAL] should never be 0.0
		}
	}
*/
// If no changes, do nothing.
	if(changeLC<0.00001 && change_crop<0.00001)
		return;
	else if(fabs(transferred_fraction-receiving_fraction)>0.0001 || fabs(change_stand-receiving_fraction)>0.0001)
			fail("Transferred landcover fractions not balanced !\n");
////////////////////////////////////////////////////////////

#define cropLUchangeCtransfer
#if defined cropLUchangeCtransfer

	double *transfer_litter_leaf, *transfer_litter_sap, *transfer_litter_heart, *transfer_litter_root, *transfer_litter_repr, *transfer_harvested_products_slow;

	transfer_litter_leaf=transfer_litter_sap=transfer_litter_heart=transfer_litter_root=transfer_litter_repr=transfer_harvested_products_slow=NULL;

	transfer_litter_leaf=new double[npft];
	transfer_litter_sap=new double[npft];
	transfer_litter_heart=new double[npft];
	transfer_litter_root=new double[npft];
	transfer_litter_repr=new double[npft];

	transfer_harvested_products_slow=new double[npft];

	double transfer_acflux_harvest=0.0;

	double transfer_cpool_fast=0.0;
	double transfer_cpool_slow=0.0;
	double transfer_wcont[NSOILLAYER]={0.0};
	double transfer_decomp_litter_mean=0.0;
	double transfer_k_soilfast_mean=0.0;
	double transfer_k_soilslow_mean=0.0;

	memset(transfer_litter_leaf,0,sizeof(double)*npft);
	memset(transfer_litter_sap,0,sizeof(double)*npft);
	memset(transfer_litter_heart,0,sizeof(double)*npft);
	memset(transfer_litter_root,0,sizeof(double)*npft);
	memset(transfer_litter_repr,0,sizeof(double)*npft);
	memset(transfer_harvested_products_slow,0,sizeof(double)*npft);

//Keep track of carbon and water in lost areas.
	gridcell.firstobj();
	while (gridcell.isobj) //Loop through stands:
	{
		double scale;

		Stand& stand=gridcell.getobj();

//		if(stand.landcover!=CROPLAND && landcoverfrac_change[stand.landcover]<0.0 || stand.landcover==CROPLAND && cropstand_change[stand.cftid]<0.0)
		if(landcoverfrac_change[stand.landcover]<0.0)
		{
//			if(stand.landcover!=CROPLAND)														
			{
				scale=-landcoverfrac_change[stand.landcover]/receiving_fraction/(double)stand.nobj;
			}
/*			else if(stand.landcover==CROPLAND)
			{
				scale=-cropstand_change[stand.cftid]/receiving_fraction/(double)stand.nobj;
			}	
*/
			stand.firstobj();
			while(stand.isobj) //Loop through Patches
			{
				Patch& patch=stand.getobj();

				Vegetation& vegetation=patch.vegetation;
				vegetation.firstobj();
				while(vegetation.isobj)
				{
					double cmass_leaf_cp=0.0, cmass_root_cp=0.0, cmass_sap_cp=0.0, cmass_heart_cp=0.0, cmass_debt_cp=0.0, cmass_ho_cp=0.0, cmass_agpool_cp=0.0, cmass_plant_cp=0.0;//bugfix 101103
					double litter_leaf_cp, litter_root_cp, litter_sap_cp, litter_heart_cp, litter_repr_cp;
					double acflux_harvest_cp = 0;
					double harvested_products_slow_cp;

					Individual& indiv=vegetation.getobj();
					Patchpft& patchpft=patch.pft[indiv.pft.id];

					cmass_leaf_cp=indiv.cmass_leaf;
					cmass_root_cp=indiv.cmass_root;
					cmass_sap_cp=indiv.cmass_sap;
					cmass_heart_cp=indiv.cmass_heart;
					cmass_debt_cp=indiv.cmass_debt;

/*					if(indiv.pft.landcover==CROPLAND)
					{
						cmass_ho_cp=indiv.cropindiv->cmass_ho;
						cmass_agpool_cp=indiv.cropindiv->cmass_agpool;
						cmass_plant_cp=indiv.cropindiv->cmass_plant;
					}
*/
					litter_leaf_cp=patchpft.litter_leaf;
					litter_root_cp=patchpft.litter_root;
					litter_sap_cp=patchpft.litter_sap;
					litter_heart_cp=patchpft.litter_heart;
					litter_repr_cp=patchpft.litter_repr;

					harvested_products_slow_cp=patch.pft[indiv.pft.id].harvested_products_slow;

//Harvest of transferred areas:
/*					if(indiv.pft.landcover==CROPLAND)
						harvest_crop(cmass_plant_cp,cmass_leaf_cp,cmass_root_cp,cmass_ho_cp,cmass_agpool_cp,
						litter_leaf_cp,litter_root_cp,acflux_harvest_cp,harvested_products_slow_cp,indiv);
					else if(patch.stand.landcover!=CROPLAND)												
*/						harvest_natural(cmass_leaf_cp,cmass_root_cp,cmass_sap_cp,cmass_heart_cp,cmass_debt_cp,	//kolla vad som händer här, både för träd och gräs !
						litter_leaf_cp,litter_root_cp,litter_sap_cp,litter_heart_cp,acflux_harvest_cp,harvested_products_slow_cp,indiv);

					gridcell.LC_updated=true;

					//In case any vegetation carbon left: (eg. cmass_root for CC3G/CC4G)
					if((cmass_leaf_cp+cmass_root_cp+cmass_sap_cp+cmass_heart_cp-cmass_debt_cp+cmass_ho_cp)!=0.0)
					{
						litter_leaf_cp+=cmass_leaf_cp;
						litter_root_cp+=cmass_root_cp;
						litter_sap_cp+=cmass_sap_cp;
						litter_heart_cp+=cmass_heart_cp-cmass_debt_cp;
/*
						if(indiv.pft.aboveground_ho)
							litter_leaf_cp+=cmass_ho_cp;
						else
							litter_root_cp+=cmass_ho_cp;
*/					}

					transfer_litter_leaf[indiv.pft.id]+=litter_leaf_cp*scale;
					transfer_litter_root[indiv.pft.id]+=litter_root_cp*scale;
					transfer_litter_sap[indiv.pft.id]+=litter_sap_cp*scale;
					transfer_litter_heart[indiv.pft.id]+=litter_heart_cp*scale;
					transfer_litter_repr[indiv.pft.id]+=litter_repr_cp*scale;

					transfer_acflux_harvest+=acflux_harvest_cp*scale;

					if(ifslowharvestpool)
						transfer_harvested_products_slow[indiv.pft.id]+=harvested_products_slow_cp*scale;

					vegetation.nextobj();
				}

//sum litter C:
				transfer_cpool_fast+=patch.soil.cpool_fast*scale;
				transfer_cpool_slow+=patch.soil.cpool_slow*scale;

//sum wcont:
				for(i=0;i<NSOILLAYER;i++)
				{
					transfer_wcont[i]+=patch.soil.wcont[i]*scale;
				}

				transfer_decomp_litter_mean+=patch.soil.decomp_litter_mean*scale;
				transfer_k_soilfast_mean+=patch.soil.k_soilfast_mean*scale;
				transfer_k_soilslow_mean+=patch.soil.k_soilslow_mean*scale;

				stand.nextobj();
			}
		}
		gridcell.nextobj();
	}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Create and kill stands:

// landcover dynamics (from updated landcoverfrac):	
	if(!lcfrac_fixed && changeLC>0.0)
	{
		for(int i=0;i<NLANDCOVERTYPES;i++)	//For all landcover types without subclasses
		{
//			if(i!=CROPLAND)
			{
				if(run[i])
				{
					if(gridcell.landcoverfrac_old[i]==0.0 && gridcell.landcoverfrac[i]>0.0)
					{
						landcover=(landcovertype)i;
						Stand& stand=gridcell.createobj(gridcell,landcover);

						pftlist.firstobj();
						while (pftlist.isobj) 
						{
							Pft& pft=pftlist.getobj();
							if(pft.landcover==i)
							{
								stand.pft[pft.id].active=true;
							}
							pftlist.nextobj();
						}
					}
					else if(gridcell.landcoverfrac_old[i]>0.0 && gridcell.landcoverfrac[i]==0.0)
					{
						gridcell.firstobj();
						while (gridcell.isobj) //Loop through stands:
						{
							Stand& stand=gridcell.getobj();
							if(stand.landcover==i)
							{
								gridcell.killobj();
							}
							else
								gridcell.nextobj();
						}
					}
				}
			}
		}
	}

//Crop stand dynamics to be put here.

//update C-pools for receiving stands:
	gridcell.firstobj();
	while (gridcell.isobj) //Loop through stands:
	{
		Stand& stand=gridcell.getobj();

//		if(stand.landcover!=CROPLAND && landcoverfrac_change[stand.landcover]>0.0 || stand.landcover==CROPLAND && cropstand_change[stand.cftid]>0.0)
		if(landcoverfrac_change[stand.landcover]>0.0)
		{
			double old_frac, added_frac, new_frac;

//			if(stand.landcover!=CROPLAND)					
			{
				old_frac=gridcell.landcoverfrac_old[stand.landcover];
				added_frac=landcoverfrac_change[stand.landcover];
				new_frac=gridcell.landcoverfrac[stand.landcover];
			}
/*			else if(stand.landcover==CROPLAND)
			{
				old_frac=gridcell.landcoverfrac_old[CROPLAND]*gridcell.cftfrac_old[stand.cftid];
				added_frac=cropstand_change[stand.cftid];
				new_frac=gridcell.landcoverfrac[CROPLAND]*gridcell.cftfrac[stand.cftid];
			}	
*/
#ifdef cropLUchangeCtransfer
			stand.firstobj();
			while(stand.isobj) //Loop through Patches
			{
				Patch& patch=stand.getobj();
//add litter C:
				for (i=0;i<npft;i++) 
				{
					Patchpft& patchpft=patch.pft[i];

					patchpft.litter_leaf=(patchpft.litter_leaf*old_frac+transfer_litter_leaf[i]*added_frac)/new_frac;
					patchpft.litter_sap=(patchpft.litter_sap*old_frac+transfer_litter_sap[i]*added_frac)/new_frac;
					patchpft.litter_heart=(patchpft.litter_heart*old_frac+transfer_litter_heart[i]*added_frac)/new_frac;
					patchpft.litter_root=(patchpft.litter_root*old_frac+transfer_litter_root[i]*added_frac)/new_frac;
					patchpft.litter_repr=(patchpft.litter_repr*old_frac+transfer_litter_repr[i]*added_frac)/new_frac;

					if(ifslowharvestpool)
						patchpft.harvested_products_slow=(patchpft.harvested_products_slow*old_frac+transfer_harvested_products_slow[i]*added_frac)/new_frac;
				}

//add soil C:
				patch.soil.cpool_fast=(patch.soil.cpool_fast*old_frac+transfer_cpool_fast*added_frac)/new_frac;
				patch.soil.cpool_slow=(patch.soil.cpool_slow*old_frac+transfer_cpool_slow*added_frac)/new_frac;

//other soil stuff:
				for(i=0;i<NSOILLAYER;i++)
					patch.soil.wcont[i]=(patch.soil.wcont[i]*old_frac+transfer_wcont[i]*added_frac)/new_frac;

				patch.soil.decomp_litter_mean=(patch.soil.decomp_litter_mean*old_frac+transfer_decomp_litter_mean*added_frac)/new_frac;
				patch.soil.k_soilfast_mean=(patch.soil.k_soilfast_mean*old_frac+transfer_k_soilfast_mean*added_frac)/new_frac;
				patch.soil.k_soilslow_mean=(patch.soil.k_soilslow_mean*old_frac+transfer_k_soilslow_mean*added_frac)/new_frac;
//add fluxes:
				patch.fluxes.report_flux(Fluxes::HARVESTC, transfer_acflux_harvest*added_frac/new_frac);

				stand.nextobj();
			}
#endif
		}
		gridcell.nextobj();
	}
#if defined cropLUchangeCtransfer
	if(transfer_litter_leaf) delete[] transfer_litter_leaf;
	if(transfer_litter_sap) delete[] transfer_litter_sap;
	if(transfer_litter_heart) delete[] transfer_litter_heart;
	if(transfer_litter_root) delete[] transfer_litter_root;
	if(transfer_litter_repr) delete[] transfer_litter_repr;
	if(transfer_harvested_products_slow) delete[] transfer_harvested_products_slow;
#endif
}
