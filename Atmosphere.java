//*****************************************************************************************
// Temperature/Atmosphere.java
// author: Non-Euclidean Dreamer
// stars with temperature, mass and radius as well as rotation information axis & daylength
//*****************************************************************************************


import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;

import javax.imageio.ImageIO;

public class Atmosphere 
{
	Planet pl;
	ArrayList<Clan> clans=new ArrayList<Clan>();
	ArrayList<Plant>plants=new ArrayList<Plant>();
	ArrayList<Animal>animals=new ArrayList<Animal>();
	static String name="try", type="png";
	static DecimalFormat df=new DecimalFormat("0000");
	static double sigma=5.67/Math.pow(10, 8),	//HelmholtzBolzmann constant
			delt=3;
	double r=600;
	double c=7000;//=specific heatcapacity*density
	double[][]temperature, population,nutrition, vegetation;
	double[][][]v;
	static double[]zro= {0,0,0};
	static Random rand=new Random();
	static boolean comment=false;
	//boolean diploid=true;
	
	//************
	//Constructors
	//************
	public Atmosphere(double radius) 
	{
		r=radius;
		temperature=new double[(int) (2*Math.PI*r)][(int) (Math.PI*r)];
		v=new double[temperature.length][temperature[0].length][2];//v[0] is angle, v[1] is norm
	}
	
	//uniform temperature t
	public Atmosphere(Planet planet,double radius,double t) 
	{
		pl=planet;
		r=radius;
		temperature=new double[(int) (2*Math.PI*r)][(int) (Math.PI*r)];
		population=new double[(int) (2*Math.PI*r)][(int) (Math.PI*r)];
		vegetation=new double[(int) (2*Math.PI*r)][(int) (Math.PI*r)];
		nutrition=new double[(int) (2*Math.PI*r)][(int) (Math.PI*r)];
		v=new double[temperature.length][temperature[0].length][2];

		for(int i=0;i<temperature.length;i++)
			for(int j=0;j<temperature[i].length;j++)
			{
				temperature[i][j]=t;
				nutrition[i][j]=20;
			}
	}

	//***************
	// Update methods
	//***************
	public void update(double phi, double psi, double in)//in=Tstar*sqrt(rstar/distance) 
	{
		//energy output E=sigmaT⁴*area
		//delT=q/c
		for(int i=0;i<temperature.length;i++)
		{	for(int j=0;j<temperature[i].length;j++) 
			{
				double dist=spherdist(phi,psi,i/r,j/r);
				if(dist<Math.PI/2)//does the spot get sun?
					temperature[i][j]+=sigma/c*(in*in*in*in*Math.cos(dist)-Math.pow(temperature[i][j], 4))*delt;//System.out.print("+");}
				else
					temperature[i][j]-=sigma/c*Math.pow(temperature[i][j], 4)*delt;//System.out.print("-");}
			}
		}		
		//A bit diffusion..bit.
		double[][]temp=temperature.clone();
		for(int i=0;i<temperature.length;i++)
			for(int j=0;j<temperature[i].length;j++) 
			{
				int jm=j-1,jp=j+1,i1=i,i2=i;
				double psired=Math.acos(Math.cos(j/r)*Math.cos(1/r)),phidiff=Math.atan2(Math.sin(1/r),Math.sin(j/r)*Math.cos(1/r));
				double[]loc0= {(i+temp.length-r*phidiff)%temp.length,r*psired},loc1= {(i+r*phidiff)%temp.length,r*psired};
				if(jm==-1) {jm=0;i1=(i+temp.length/2)%temp.length;}
				if(jp==temp[0].length) {jp--;i2=(i+temp.length/2)%temp.length;}
				temperature[i][j]=0.25*(temp[i1][jm]+temp[i2][jp]+FlowField.average(temp, loc1)+FlowField.average(temp, loc0));
			}
				
		
	}
	
	public void update(ArrayList<Star> stars,double t) 
	{
		for(int i=0;i<temperature.length;i++)
		{	for(int j=0;j<temperature[i].length;j++) 
			{
			temperature[i][j]-=sigma/c*Math.pow(temperature[i][j], 4)*t;//}
			}
		}
		for(Star st:stars)
		{
			double[] ph=phipsi(st); //where on the map is the star in zenit?
			double d=norm(subtract(st.loc,pl.loc));
			for(int i=0;i<temperature.length;i++)
			{	for(int j=0;j<temperature[i].length;j++) 
				{
				double dist=spherdist(ph[0]+pl.greenwich,ph[1],i/r,j/r);
				
					if(dist<Math.PI/2)
						{
							double in=st.temperature*Math.sqrt(st.radius/d);
							temperature[i][j]+=sigma/c*(in*in*in*in*Math.cos(dist))*t;
						}
				}
			}
		}
		
		
		double[][]temp=temperature.clone();
		for(int i=0;i<temperature.length;i++)
			for(int j=0;j<temperature[i].length;j++) 
			{
				int jm=j-1,jp=j+1,i1=i,i2=i;
				double psired=Math.acos(Math.cos(j/r)*Math.cos(1/r)),phidiff=Math.atan2(Math.sin(1/r),Math.sin(j/r)*Math.cos(1/r));
				double[]loc0= {(i+temp.length-r*phidiff)%temp.length,r*psired},loc1= {(i+r*phidiff)%temp.length,r*psired};
				if(jm==-1) {jm=0;i1=(i+temp.length/2)%temp.length;}
				if(jp==temp[0].length) {jp--;i2=(i+temp.length/2)%temp.length;}
				temperature[i][j]=0.25*(temp[i1][jm]+temp[i2][jp]+FlowField.average(temp, loc1)+FlowField.average(temp, loc0));
				//Chances for spawning a clan
				//if(population[i][j]<0.01)	if(rand.nextDouble()<0.000000000005*Math.sin(j/r)*Math.pow(40-Math.abs(temperature[i][j]-295),3))clans.add(new Clan(this, new double[] {i,j}));
				//Chances for spawning a plant
			//	if(diploid)
				{
					if(vegetation[i][j]+population[i][j]<0.01)	if(rand.nextDouble()<0.00000000024*Math.sin(j/r)*Math.pow(40-Math.abs(temperature[i][j]-295),3)/Math.log(Math.E+plants.size()+animals.size()))plants.add(new Plant(this, new double[] {i,j},1));

				}
			//	else 
				{
				if(vegetation[i][j]<0.01)	if(rand.nextDouble()<0.00000000008*Math.sin(j/r)*Math.pow(40-Math.abs(temperature[i][j]-295),3)/Math.log(Math.E+plants.size()))plants.add(new Plant(this, new double[] {i,j}));
				//Chances for spawning a clan
				if(population[i][j]<0.01)	if(rand.nextDouble()<0.0000000003*Math.sin(j/r)*Math.pow(40-Math.abs(temperature[i][j]-295),3)/Math.log(Math.E+animals.size()))animals.add(new Animal(this, new double[] {i,j}));
			}}
			
	/*	if(diploid)
		{
			if(plants.size()+animals.size()==0)plants.add(new Plant(this,new double[] {rand.nextDouble()*temperature.length,rand.nextDouble()*temperature[0].length},1));
		}
		else {
		if(plants.size()==0)plants.add(new Plant(this,new double[] {rand.nextDouble()*temperature.length,rand.nextDouble()*temperature[0].length}));
		if(animals.size()==0)animals.add(new Animal(this,new double[] {rand.nextDouble()*temperature.length,rand.nextDouble()*temperature[0].length}));}*/
//plants
		for(int i=plants.size()-1;i>-1;i--)
		{
			Plant pl=plants.get(i);
			int x=(int)pl.loc[0],y=(int)pl.loc[1];
		//	System.out.println("energy="+pl.energy+", minres="+pl.minreservoir+", size="+pl.size+", maxsize="+pl.maxsize+", leaves="+pl.leaves);
			double tmp=temperature[x][y],  diff=Math.abs(tmp-pl.idealtemp);
			int sg=1;
			if(tmp<pl.idealtemp)sg=0;
			double win=(tmp*pl.leavefficiency*Math.min(1,Math.pow((1+pl.size)/(1+vegetation[x][y]-pl.size),5.2)/1000)-diff/Math.pow(1+pl.durability[2+sg],2)/6000);
			if(pl.leaves>0)
				{
				  win*=pl.leaves;
				  if(win>0)
					 {
				 		pl.energy+=win;
				  		if(rand.nextDouble()<pl.beh[0])pl.growLeaf();
				  }
				  else if(rand.nextDouble()<pl.beh[1]) {pl.leaves--;//System.out.println("losing leaves");
				  nutrition[x][y]+=(1+pl.leavefficiency)*0.15;}
				  else pl.energy+=win;
				}
			else
			{
				if(win>0&&rand.nextDouble()<pl.beh[3])pl.growLeaf();
			}
			if(pl.size>0)
			pl.energy-=pl.size*(diff*diff)/Math.pow(1+pl.durability[sg],2)/400;
			else
				pl.energy-=diff/Math.pow(1+pl.durability[sg+4], 2)/3000;
			if( pl.energy<0) 
			{
				nutrition[x][y]+=pl.size*pl.growenergy();
				vegetation[(int) pl.loc[0]][(int) pl.loc[1]]-=pl.size;
				plants.remove(i);
			}
			else if(pl.energy>pl.minreservoir)
			{
				if(rand.nextDouble()<pl.beh[4]/delt)
				{
					if(pl.size<pl.maxsize)pl.grow();
					else 
					{
					//	System.out.println("plant fullsized and ready");
						if(pl.diploid)animals.add(pl.baby());
						else plants.add(pl.seed());
					}
				}	
			}
		}
		//animals they have two moves per round. eat&then move or multiply or grow
		for(int i=animals.size()-1;i>-1;i--)
		{
			Animal an=animals.get(i);
			double oldenergy=an.energy;
			double tmp=temperature[(int)an.loc[0]][(int)an.loc[1]],diff=Math.abs(tmp-an.idealtemp);
			int sg=1;
			if(tmp<an.idealtemp)sg=0;
			int x=(int)an.loc[0],y=(int)an.loc[1];
			
			if(vegetation[x][y]*(an.leavestomach+an.plantstomach)>nutrition[x][y]*an.carrionstomach)
			{
				boolean found=false;
				int j=0;
				while(!found&j<plants.size())
				{
					Plant pl=plants.get(j);
					if (x==(int)pl.loc[0]&&y==(int)pl.loc[1])
					{
						double r=rand.nextDouble();
						//if(pl.leaves+pl.size==0)an.eatSeed(pl);
						
						 if(r<an.leavestomach/(1+an.plantstomach+an.leavestomach)&&pl.leaves>0)an.eatLeaf(pl);
						else if(r<(an.leavestomach+an.plantstomach)/(1+an.plantstomach+an.leavestomach)&&pl.size>0)an.eatPlant(pl);
						else		an.eatBattery(pl);
						found=true;
						
					}
					else j++;
				}
				
			}
			else an.eatNutrition();
			
			an.energy-=an.livenergy(diff,sg);

			if(an.energy<0) {System.out.println("demise of animal("+an.loc[0]+","+an.loc[1]+")");population[(int)an.loc[0]][(int)an.loc[1]]-=an.size*an.growenergy();animals.remove(i);nutrition[x][y]+=an.size/2;}
			else if(population[(int)an.loc[0]][(int)an.loc[1]]>an.size) {if(an.energy<oldenergy)
			{
				double r=rand.nextDouble();
				if(r<an.beh[0])an.moveTemp();
				else if(r<an.beh[0]+an.beh[1]) an.movetoFood();
			}
			else if(an.energy>an.minreservoir)
			{
				double r=rand.nextDouble();
				if(r<an.beh[2]) {
				if(an.size>an.maxsize) {if(an.diploid)plants.add(an.seedout());else animals.add(an.seed());}
				else an.grow();}
			//	else if(r<0.6)an.growstomach();
			}}
			else
			{
				if(an.energy>an.minreservoir)
				{
					double r=rand.nextDouble();
					if(r<an.beh[2]) {
					if(an.size>an.maxsize) {if(an.diploid)plants.add(an.seedout());else animals.add(an.seed());}
					else an.grow();}
				//	else if(r<0.6)an.growstomach();
				}
				else if(an.energy<oldenergy)
				{
					double r=rand.nextDouble();
					if(r<an.beh[0])an.moveTemp();
					else if(r<an.beh[0]+an.beh[1]) an.movetoFood();
				}
			}

		}
		
			//clans
			for(int i=clans.size()-1;i>-1;i--)
			{
				Clan cl=clans.get(i);
				double tmp=temperature[(int)cl.loc[0]][(int)cl.loc[1]];
				if(tmp>cl.livzone[1]+cl.comfortzone[1]||tmp<cl.livzone[0]+cl.comfortzone[0]||population[(int)cl.loc[0]][(int)cl.loc[1]]>Math.sin(cl.loc[1]/r)*100*cl.size) {System.out.println("demise of clan("+cl.loc[0]+","+cl.loc[1]+")");population[(int)cl.loc[0]][(int)cl.loc[1]]-=cl.size;clans.remove(i);}
				else if(tmp<cl.comfortzone[1]&&tmp>cl.comfortzone[0]&&population[(int)cl.loc[0]][(int)cl.loc[1]]<Math.sin(cl.loc[1]/r)*10*cl.size)
				{
					if(rand.nextDouble()<.25&&rand.nextDouble()*cl.size>3.5)clans.add(cl.split());
					else if(rand.nextDouble()<0.7) cl.evolve(rand.nextDouble()*0.8);
				}
				else
				{
					cl.move();
				}
			}
	
	
	//decompose
	for(int i=0;i<nutrition.length;i++)
		for(int j=0;j<nutrition[i].length;j++)nutrition[i][j]/=Math.pow(temperature[i][j]/375,16)+1;
	}
	//for wind
	private void advect(double t) 
	{
		int x=temperature.length,y=x/2;
		double[][][]vc=v.clone();
		for(int i=0;i<x;i++)
		{	for(int j=0;j<y;j++) 
			if(v[i][j][1]>0.01)
			{
				double psi=(j+0.1)/r,vn=vc[i][j][1]/r,alpha=vc[i][j][0],psinew=Math.acos(Math.cos(psi)*Math.cos(vn)+Math.sin(psi)*Math.sin(vn)*Math.cos(alpha)) ,
						delphi=Math.acos((Math.cos(vn)-Math.cos(psinew)*Math.cos(psi))/(Math.sin(psinew)*Math.sin(psi))), alphanew=Math.acos((Math.cos(vn)*Math.cos(psinew)-Math.cos(psi))/(Math.sin(vn)*Math.sin(psinew))),
						inew=(i+r*delphi+x)%x, jnew=Math.min(Math.max(0.001, r*psinew),y-0.1);
				double[]lc= {inew,jnew};
				
				v[i][j]=FlowField.average(v,lc);
				v[i][j][0]+=alpha-alphanew;
				if(Double.isNaN(v[i][j][1])||Double.isInfinite(v[i][j][1])) {v[i][j][1]=0;v[i][j][0]=0;}//System.out.println("advect error");
				else temperature[i][j]=FlowField.average(temperature, lc);
			
			}
		}
	}

	public void updateWind(ArrayList<Star> stars,double t) 
	{
		for(int i=0;i<temperature.length;i++)
		{	for(int j=0;j<temperature[i].length;j++) 
			{
			temperature[i][j]-=sigma/c*Math.pow(temperature[i][j], 4)*t;//}
			}
		}
		if(comment) {System.out.println("Cool Down");	FlowField.print(temperature);}
		for(Star st:stars)
		{
			double[] ph=phipsi(st); //where on the map is the star in zenit?
			double d=norm(subtract(st.loc,pl.loc));
			for(int i=0;i<temperature.length;i++)
			{	for(int j=0;j<temperature[i].length;j++) 
				{
				double dist=spherdist(ph[0]+pl.greenwich,ph[1],i/r,j/r);
				
					if(dist<Math.PI/2)
						{
							double in=st.temperature*Math.sqrt(st.radius/d);
							temperature[i][j]+=sigma/c*(in*in*in*in*Math.cos(dist))*t;
						}
					
					//Chances for spawning a clan
					//	if(population[i][j]<0.01)	if(rand.nextDouble()<0.000000000005*Math.sin(j/r)*Math.pow(40-Math.abs(temperature[i][j]-295),3))clans.add(new Clan(this, new double[] {i,j}));
					}}
				/*
					//clans
					for(int i=clans.size()-1;i>-1;i--)
					{
						Clan cl=clans.get(i);
						double tmp=temperature[(int)cl.loc[0]][(int)cl.loc[1]];
						if(tmp>cl.livzone[1]+cl.comfortzone[1]||tmp<cl.livzone[0]+cl.comfortzone[0]||population[(int)cl.loc[0]][(int)cl.loc[1]]>Math.sin(cl.loc[1]/r)*100*cl.size) {System.out.println("demise of clan("+cl.loc[0]+","+cl.loc[1]+")");population[(int)cl.loc[0]][(int)cl.loc[1]]-=cl.size;clans.remove(i);}
						else if(tmp<cl.comfortzone[1]&&tmp>cl.comfortzone[0]&&population[(int)cl.loc[0]][(int)cl.loc[1]]<Math.sin(cl.loc[1]/r)*10*cl.size)
						{
							if(rand.nextDouble()<.25&&rand.nextDouble()*cl.size>3.5)clans.add(cl.split());
							else if(rand.nextDouble()<0.7) cl.evolve(rand.nextDouble()*0.8);
						}
						else
						{
							cl.move();
						}
					}
			*/
		}
		if(comment) {System.out.println("Warm up");	FlowField.print(temperature);}

		//advect
				advect(delt);
				
				if(comment) {	System.out.println("advect");FlowField.print(temperature);FlowField.print(v);}
				corriolis(0.5);
				if(comment) {System.out.println("Wind up");	FlowField.print(v);}
				relax(pl.radius/r);
				if(comment) {System.out.println("relax");	FlowField.print(v);}
		double[][]temp=temperature.clone();
		
		
		for(int i=0;i<temperature.length;i++)
			for(int j=0;j<temperature[i].length;j++) 
			{
				int jm=j-1,jp=j+1,i1=i,i2=i;
				double psired=Math.acos(Math.cos(j/r)*Math.cos(1/r)),phidiff=Math.atan2(Math.sin(1/r),Math.sin(j/r)*Math.cos(1/r));
				double[]loc0= {(i+temp.length-r*phidiff)%temp.length,r*psired},loc1= {(i+r*phidiff)%temp.length,r*psired};
				if(jm==-1) {jm=0;i1=(i+temp.length/2)%temp.length;}
				if(jp==temp[0].length) {jp--;i2=(i+temp.length/2)%temp.length;}
				temperature[i][j]=0.25*(temp[i1][jm]+temp[i2][jp]+FlowField.average(temp, loc1)+FlowField.average(temp, loc0));
				//Chances for spawning a clan
			//	if(temperature[i][j]<313&&temperature[i][j]>263)
				//	if(rand.nextDouble()<0.000000000005*Math.sin(j/r)*Math.pow(40-Math.abs(temperature[i][j]-295),3))clans.add(new Clan(this, new double[] {i,j}));
			}
	if(comment)	{	System.out.println("diffuse");	FlowField.print(temperature);FlowField.print(v);}
	}

	private void relax(double delx)
	{
		double max=10,bound=0.001;
		int l=0,k=0;
		//System.out.println("i="+f.vectorfield[0][0].length);
		int x=temperature.length,y=x/2;

		double[][][]q=new double[2][x][y];
		double[][]div=divergence();
		while(max>bound&&l<100)//500l<400&&
		{
			
			//poisson:
		
			max=0;
			for(int i=0;i<x;i++)
				for(int j=0;j<y;j++)
				{
					int jm=j-1,jp=j+1,i1=i,i2=i;
					double psi=(j+0.1)/r,psired=Math.acos(Math.cos(psi)*Math.cos(1/r)),phidiff=Math.atan2(Math.sin(1/r),Math.sin(psi)*Math.cos(1/r));
					double[]loc0= {(i+x-r*phidiff)%x,Math.max(y-.0001, r*psired)},loc1= {(i+r*phidiff)%x,Math.max(y-.0001, r*psired)};
					if(jm==-1) {jm=0;i1=(i+x/2)%x;}
					if(jp==y) {jp--;i2=(i+x/2)%x;}
					q[1-k][i][j]=0.25*(q[k][i1][jm]+q[k][i2][jp]+FlowField.average(q[k], loc1)+FlowField.average(q[k], loc0)-0.01*div[i][j]);
					max=Math.max(max, Math.abs(q[0][i][j]-q[1][i][j]));
				}
			l++;k=l%2;
		}
		//		System.out.println("q");
		//	FlowField.print(q[k]);
		for(int i=0;i<x;i++)
			for(int j=0;j<y;j++)
				{
					double vx=Math.sin(v[i][j][0])*v[i][j][1],vy=Math.cos(v[i][j][0])*v[i][j][1];
					int jm=j-1,jp=j+1,i1=i,i2=i;
					if(jm==-1) {jm=0;i1=(i+x/2)%x;}
					if(jp==y) {jp--;i2=(i+x/2)%x;}
					vy-=(q[k][i2][jp]-q[k][i1][jm])/2;
					
					i1=(i+x-1)%x; i2=(i+1)%x;
					vx-=(q[k][i2][j]-q[k][i1][j])/2/Math.max(0.01, Math.sin(j/r));
					v[i][j][0]=Math.atan2(vx, vy);
					v[i][j][1]=Math.sqrt(vx*vx+vy*vy);
					if(Double.isNaN(v[i][j][1])||Double.isInfinite(v[i][j][1]))System.out.println("relax error");

				}
		
	}

	private double[][] divergence() 
	{
		int x=temperature.length,y=x/2;
		double delx=pl.radius/r;
	double [][]out=new double[x][y];//	System.out.println("div=");FlowField.print(v);
		for(int i=0;i<x;i++)
			for(int j=0;j<y;j++)
			{
				int jm=j-1,jp=j+1,i1=i,i2=i;
				double psired=Math.acos(Math.cos(j/r)*Math.cos(1/r)),phidiff=Math.atan2(Math.sin(1/r),Math.sin(j/r)*Math.cos(1/r));
				double[]loc0= {(i+x-r*phidiff)%x,r*psired},loc1= {(i+r*phidiff)%x,r*psired};
				if(jm==-1) {jm=0;i1=(i+x/2)%x;}
				if(jp==y) {jp--;i2=(i+x/2)%x;}
				double[] vp=FlowField.average(v, loc1), vm=FlowField.average(v, loc0);
				out[i][j]=0.5*(-v[i1][jm][1]*Math.cos(v[i1][jm][0])+v[i2][jp][1]*Math.cos(v[i2][jp][1])+vp[1]*Math.sin(vp[0])-vm[1]*Math.sin(vm[0]));	
				
				//FlowField.print(loc1);
			}
			
		return out;
	}

	//distance between spots on sphere
	private double spherdist(double phi, double psi, double phi1, double psi1) 
	{
		return Math.acos(Math.sin(psi)*Math.sin(psi1)*Math.cos(phi-phi1)+Math.cos(psi)*Math.cos(psi1));
	}
	
	//where on the sphere is the star in zenit
	private double[] phipsi(Star st) 
	{
		double[]dir=subtract(st.loc,pl.loc);
		
		double norm=norm(dir),cospsi=dot(pl.axis,dir)/norm;//System.out.println(pl.phinull[0]+","+pl.phinull[1]+","+pl.phinull[2]);//Math.signum(dot(pl.phi90,dir)));
		return new double[] {Math.signum(dot(pl.phi90,dir))*Math.acos(dot(pl.phinull,dir)/norm)/cospsi ,Math.acos(cospsi)};
	}


	//*******************
	//Draw the atmosphere
	//*******************
	public void print(String name,String type)
	{
		System.out.println(temperature[10][100]);
		File file=new File(name);
		BufferedImage image=new BufferedImage(temperature.length,temperature[0].length,BufferedImage.TYPE_4BYTE_ABGR);
		for(int i=0;i<image.getWidth();i++)
			for(int j=0;j<image.getHeight();j++)
				image.setRGB(i,j,colorcode(i,j));
		
		try {
				ImageIO.write(image, type, file);
			}	catch (IOException e) {	System.out.println("IOException: Problems saving file "+name);	e.printStackTrace();}
	}

	int colorcode(int i,int j) 
	{
		double d=temperature[i][j];
		int step=25,start=290-4*step;
		int[][]cl= {{128,128,128},{255,255,255},{128,128,255},{0,0,255},{0,255,255},{0,255,0},{255,255,0},{255,128,0},{255,0,0},{0,0,0}};
		for(int k=0;k<cl.length-1;k++)
		{
			if(d<start)
			{
				double st=(start-d)/step;
				try{return new Color((int)(cl[k][0]*st+cl[k+1][0]*(1-st)),(int)(cl[k][1]*st+cl[k+1][1]*(1-st)),(int) (cl[k][2]*st+cl[k+1][2]*(1-st))).getRGB();}
				catch(IllegalArgumentException e) {return Color.black.getRGB();}
			}
			start+=step;
		}
		int n=cl.length-1;
		return new Color(cl[n][0],cl[n][1],cl[n][2]).getRGB();
	}

	
	//**************
	// vector methods
	//**************
	static double norm(double[] v) 
	{
		return Math.sqrt(dot(v,v));
	}

	private static double dot(double[] v, double[] w) 
	{
		double out=0;
		for(int i=0;i<v.length;i++)out+=v[i]*w[i];
		return out;
	}

	private double[] subtract(double[] v, double[] w)
	{
		double[]out=new double[v.length];
		for(int i=0;i<v.length;i++)
			out[i]=v[i]-w[i];
		return out;
	}
	
	public void setC(double c1)
	{
		c=c1;
	}
	
	public double[] gradient(int i,int j)
	{
		int jp=Math.min(j+1,temperature[0].length-1),jm=Math.max(j-1,0),im=i-1,ip=i+1;
		if(im<0)im=temperature.length-1;
		if(ip==temperature.length)ip=0;
		
		return new double[] {(temperature[ip][j]-temperature[im][j])/Math.sin((j+0.1)/r),temperature[i][jp]-temperature[i][jm]};
		
	}

	public void setT(double t) 
	{
		for(int i=0;i<temperature.length;i++)
			for(int j=0;j<temperature[i].length;j++)
				temperature[i][j]=t;
	}
	
	public void corriolis(double strength)
	{
		for(int i=0;i<temperature.length;i++)
			for(int j=0;j<temperature[i].length;j++)
			{
				double[] gr=gradient(i,j);
				//simulation the pressure bias, although not correctly
				double x=Math.sin(v[i][j][0])*v[i][j][1]-strength*gr[0],
						y=Math.cos(v[i][j][0])*v[i][j][1]-strength*gr[1];
				//corriolis
				double factor=2*Math.cos(j/r)/pl.sidday,
						ax=-factor*y, ay=factor*x;
				x=x+ax;y=y+ay;
				v[i][j][0]=Math.atan2(x,y);//Math.pow(Math.sin(j/r),2)*strength/pl.sidday+
				v[i][j][1]=Math.sqrt(x*x+y*y);
				if(Double.isNaN(v[i][j][1])||Double.isInfinite(v[i][j][1])){v[i][j][0]=0;v[i][j][1]=0;}//System.out.println("corriolis error");

			}
	}
	
	//Test method from before merging it with orbits
	/*public static void main(String[] args)
	{
		int t=-100;
		double gamma=Math.PI/6,s=0;
	//	Atmosphere a=new Atmosphere(new Planet(zro.clone(),zro.clone(),zro.clone(),300,1),300,273);
		
		while(true) 
		{
			System.out.println(t);
			a.update(t*0.1+Math.atan2(Math.sin(s), Math.cos(gamma)*Math.cos(s)),Math.acos(-Math.sin(gamma)*Math.cos(s)),5900/15.0);
			if(t>-1)
			a.print(name+df.format(t)+"."+type, type);
			t++;
			s+=0.003;
		}
	}*/
}
