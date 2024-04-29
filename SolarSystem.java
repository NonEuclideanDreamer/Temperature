//************************************************************************
// Temperature/SolarSystem.jave
// Main Method for simulating the orbits and temperatures in a star system

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File; 
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import javax.imageio.ImageIO;

public class SolarSystem  
{  
	ArrayList<Planet> planets;
	ArrayList<Star> stars;
	static String name="oneSun",type="png";
	static DecimalFormat df=new DecimalFormat("0000");
	static double[][]zBuffer=new double[2560][1440];//
	//physical constants in km,kg,h units
	static double g=8.6E-13;
	static int dim=3; 
	static double rulescope=dim-1;
	static double[]zro= {0,0,0};
	static double metric=2, delt=20;
	static int steps=1,black=Color.black.getRGB(), white=Color.white.getRGB(),k=1;//steps:How many iteration steps "+delt" between images	
	static double rsun=690000,msun=2E30,au=1.5E8,
			time=2;
	static boolean atm=false, big=false, doble=false;
	static Random rand=new Random();
	
	//*************
	// Constructor
	//*************
	public SolarSystem() 
	{
		planets=new ArrayList<Planet>(); 
		stars=new ArrayList<Star>();
	}
	void add(Planet planet) 
	{
		planets.add(planet);	
	}

	void add(Star star) 
	{
		stars.add(star);
	}
	
	//************
	// Main Method 
	//************
	public static void main(String[]args)
	{ 
		SolarSystem sys=starshape();//twoAndOne();//figureight();//oneSun();//trefoil();//wild();//binary();//ArmaidaPuck();//symmetric();//circleSuns(1000);//
		if (big)zBuffer=new double[2560][1440];
		int counter=0;
		Atmosphere a=sys.planets.get(0).atm; 
	if(atm) {	
		a.setT(200);
		if(doble) {
		a=sys.planets.get(1).atm; 
		a.setT(200);
		//a.diploid=false;
		}}
		BufferedImage image=new BufferedImage(1440*k,2560*k,BufferedImage.TYPE_4BYTE_ABGR);
		if(!big &&atm)		image=new BufferedImage(/*2560,1440*/a.temperature.length,a.temperature.length,BufferedImage.TYPE_4BYTE_ABGR);
		//else image=new BufferedImage(zBuffer.length,zBuffer[0].length,BufferedImage.TYPE_4BYTE_ABGR);
		for(int j=0;j<0;j++)
		{ sys.log(time);
			System.out.println(counter); 
			for(int i=0;i<steps*24;i++)	
			{
				sys.update(delt);time=time+delt/24.0; 
			}
			
		}
		//delt/=10;
		//atm=true;
		int k=0;
		while(true)//k<400)  
		{ //sys.log(time);
			k++;
			System.out.println(counter); 
			for(int i=0;i<steps;i++)	sys.update(delt);
			time+=delt*steps/24.0;
			if(big&&!doble)printwrite(image,sys,0,counter);
			else if(atm) print(image,sys,0,counter);
			else print(image,sys,counter);
			counter++; 
		//	System.out.println("temp="+a.temperature[5][5]+", v=("+a.v[5][5][0]+","+a.v[5][5][1]+")");
			
		}
	}
	
	private static SolarSystem circle() 
	{
		SolarSystem out=new SolarSystem();
		double a=au,v=1E7,mass=1E30;
		out.add(new Star(new double[] {-a,0,0},new double[] {0,v*mass,0},2E7,mass,0));
		out.add(new Star(new double[] {a,0,0},new double[] {0,-v*mass,0},2E7,mass,0));
		return out;
	}
	private static SolarSystem circleSuns(int n) 
	{ 
		double r=150*au;
		SolarSystem out=new SolarSystem();
		for(int i=0;i<n;i++)
		{
			double phi=Math.PI*2*i/Math.PI,mass=1E32;
			out.add(new Star(new double[] {r*Math.cos(phi),r*Math.sin(phi),0}, new double[] {Math.sin(phi)*mass*1E7,-Math.cos(phi)*mass*1E7,0}, 3000000,mass,0));
			r*=1.0002;
		}
		return out;
	}
	//************************
	// Specific Solar Systems
	//************************
	public static SolarSystem trefoil()//2116
	{
		SolarSystem sys=new SolarSystem();
		double mass=msun,a=1.5,v=2E4;
		sys.add(new Star(new double[] {-2*a*au,0,0},new double[] {0,-2*v*mass,0},rsun,mass,2000));
		
		sys.add(new Star(new double[] {a*au,a*Math.sqrt(3)*au,0*au},new double[] {-Math.sqrt(3)*v*mass,v*mass,0},rsun*1.5,mass,4000));
	
		sys.add(new Star(new double[]{a*au,-a*Math.sqrt(3)*au,0*au},new double[] {Math.sqrt(3)*v*mass,v*mass,0},rsun*0.75,mass,7000));
		
		mass=8E24;
		sys.add(new Planet(new double[]{2.4*au,-.7*au,-1.4*au},new double[] {1.5E4*mass,1E4*mass,1.7E4*mass},new double[] {0.1,0.4,1},6000,mass,Color.green));
		sys.add(new Planet(new double[]{2.5*au,-.8*au,-0.7*au},new double[] {1.66E4*mass,4.29E4*mass,-8.48E4*mass},new double[] {0.1,0.4,1},6000,mass,Color.cyan));

		return sys;
	}
	public static SolarSystem figureight()//2116
	{
		SolarSystem sys=new SolarSystem();
		double mass=msun,a=3*au,v=6.2E4*mass;
		sys.add(new Star(new double[] {-0.97000436*a, 0.24308753*a,0},new double[] {0.4662036850*v, 0.4323657300*v,0},rsun*1.25,mass,7000));
		
		sys.add(new Star(new double[] {0,0,0},new double[] {-0.93240737*v, -0.86473146*v,0},rsun*1.25,mass,7000));
	
		sys.add(new Star(new double[]{0.97000436*a, -0.24308753*a,0*au},new double[] {0.4662036850*v, 0.4323657300*v,0},rsun*1.25,mass,7000));
		
		mass=8E24;
		for(int i=0;i<256;i+=1)//127.3 (80,-90,115,120,120,225 also potential)
		{sys.add(new Planet(new double[]{au*2.170878,0*au,-0*au},new double[] {0E4*mass,(2-i/16.0)*1E3*mass,0E4*mass},new double[] {0.1,0.4,1},6000*Math.pow(Math.sin(i*Math.PI/256.0),0.2),mass,new Color(i,i,128)));
		sys.add(new Planet(new double[]{au*2.046,0*au,-0*au},new double[] {0E4*mass,-(48+i/32.0)*1E3*mass,0E4*mass},new double[] {0.1,0.4,1},6000*Math.pow(Math.sin(i*Math.PI/256.0),0.2),mass,new Color(i,128,i)));
		sys.add(new Planet(new double[]{0*au,1.273*au,-0*au},new double[] {0E4*mass,-(i/32.0)*1E3*mass,0E4*mass},new double[] {0.1,0.4,1},6000*Math.pow(Math.sin(i*Math.PI/256.0),0.2),mass,new Color(128,i,i)));
		}
		return sys;
	}
	public static SolarSystem twoAndOne()//2116
	{
		SolarSystem sys=new SolarSystem();
		double mass=msun,b=au/2,n=4.075,e=0,E=0,factor=Math.pow(n*n*2,1.0/3), a=b*factor,w=Math.sqrt(g*(1+e)*(1-e)/b*mass)*mass*0.5,v=Math.sqrt(g*(1+E)*(1-E)/a*mass)*mass*0.339;//w can be negative or positive
		sys.add(new Star(new double[] {b*(1-e)-a*(1-E),0,0},new double[] {0*v, v+w,0},rsun*1,mass,8000));
		
		sys.add(new Star(new double[] {-b*(1-e)-a*(1-E),0,0},new double[] {-0.*v, v-w,0},rsun*1,mass,7000));
	
		sys.add(new Star(new double[]{2*a*(1-E), -0*a,0*au},new double[] {0*v, -2*v,0},rsun*2,mass,5800));
		
		mass=8E24;
		/*for(int i=0;i<256;i+=1)
		{sys.add(new Planet(new double[]{-a,au*1.3,-0*au},new double[] {1E3*(125+(i-128)/25.0)*mass,v/msun*mass,0E4*mass},new double[] {0.1,0.1,1},6000,mass,new Color(i,128,128)));
		sys.add(new Planet(new double[]{3*au,au*0,-0*au},new double[] {1E3*(i-128)/25.0*mass,1.25E5*mass,0E4*mass},new double[] {0.1,0.4,1},6000,mass,new Color(128,i,128)));
		sys.add(new Planet(new double[]{1*au,au*0,-0*au},new double[] {0*mass,1E3*(-124+(i-128)/25.0)*mass,0E4*mass},new double[] {0.1,0.4,1},6000,mass,new Color(128,128,i)));
		}*/		
		sys.add(new Planet(new double[]{1*au,au*0,-0*au},new double[] {0*mass,1E3*(-124+(-128)/25.0)*mass,0E4*mass},new double[] {0.1,0.4,1},6000,mass,new Color(0,255,0)));
		
		return sys;
	}
	public static SolarSystem starshape()//2116
	{
		SolarSystem sys=new SolarSystem();
	//	double mass=msun,a=3*au,d=0.07,v=mass*6.4E4,w=mass*1.21E4,b=d*a,c=(d+1)*a;//n=2
		double mass=msun,a=3*au,d=0.287,v=mass*7.29E4,w=mass*0.407E4,b=d*a,c=(d+1)*a;//n=3
		sys.add(new Star(new double[] {-c,0,0},new double[] {0,-w,0},rsun*1,mass,7000));
		
		sys.add(new Star(new double[] {b,0,0},new double[] {0, v+w,0},rsun*1,mass,7000));
	
		sys.add(new Star(new double[]{a,0,0},new double[] {0, -v,0},rsun,mass,7000));
		
		mass=8E24;
	
		//sys.add(new Planet(new double[]{-1.9*au,au*0.6,-0*au},new double[] {1E3*(-60.212)*mass,-0.19999E4*mass,0E4*mass},new double[] {0.1,0.1,1},6000,mass,Color.black));
		for(int i=0;i<256;i+=1) {
		sys.add(new Planet(new double[]{au/100*(-58),au*0.6,-0*au},new double[] {1E3*(-10+(i-130)/64.0)*mass,-1.51E3*mass,0E4*mass},new double[] {0.1,0.4,1},6000*Math.pow(Math.sin(i*Math.PI/256.0),0.2),mass,new Color(i,i,128)));
		sys.add(new Planet(new double[]{3*au,au*1,-0*au},new double[] {2E3*(34.665+(i)/1024.0)*mass,(-0)*1E3*mass,0*mass},new double[] {0.1,0.4,1},6000*Math.pow(Math.sin(i*Math.PI/256.0),0.2),mass,new Color(i,128,i)));
		sys.add(new Planet(new double[]{au/100*(-29),au*0.6,-0*au},new double[] {1E3*(-10+(i-130)/64.0)*mass,(-1.01+(i-130)/256.0)*1E3*mass,0*mass},new double[] {0.1,0.4,1},6000*Math.pow(Math.sin(i*Math.PI/256.0),0.2),mass,new Color(128,i,i)));

		}
		return sys;
	}
	public static SolarSystem ArmaidaPuck()
	{
		SolarSystem sys=new SolarSystem(); 
		double 
		mass=8.31E27;
		sys.add(new Star(new double[] {-5.25E10*.592,0,0},new double[] {0*mass,-857*mass,.0*mass},13.1E6,mass,10110));//Armaida
		mass=5.27E26;
		sys.add(new Star(new double[] {8.28E11*.592,0,0},new double[] {-0*mass,13510*mass,-.00*mass},558.000,mass,12549));//Puck
		
		mass=9.1388E22;
		sys.add(new Planet(new double[]{-5.25E10*.592+1E9,1E10,0},new double[] {-2E4*mass,(-857-1.72E2)*mass,(0)*mass},new double[] {0.1,0.1,1},1887.398,mass,Color.green));//Saimun
		mass=7.1388E22;
		sys.add(new Planet(new double[]{-5.25E10*.592-1E9,1E10,0},new double[] {+2E4*mass,(-857-1.72E2)*mass,(0)*mass},new double[] {0.1,0.3,0.2},1887.398*1.5,mass,Color.cyan));//Reimun

		return sys;
	} 
	
	public static SolarSystem wild()//cooks at 3069
	{
		SolarSystem sys=new SolarSystem(); 
		double mass=msun;
		sys.add(new Star(new double[] {-2*au,0,0},new double[] {-5E3*mass,-2E4*mass,.006*mass},rsun,mass,5800));
		mass=1*msun;
		sys.add(new Star(new double[] {2*au,0,0},new double[] {-8E2*mass,2.24E4*mass,-.003*mass},1*rsun,mass,9000));
		
		mass=.5*msun;
		sys.add(new Star(new double[]{0,0,2*au},new double[] {1.62E4*mass,-.48E4*mass,.006*mass},rsun*0.5,mass,3660));
		
		mass=8E24;
		sys.add(new Planet(new double[]{-2*au,0,-2*au/*-1*/},new double[] {-1E4*mass,-1.7E4*mass,(.004)*mass},new double[] {0.1,0.1,1},6000,mass,Color.green));//4860

		return sys;
	} 
	public static SolarSystem threed()
	{
		SolarSystem sys=new SolarSystem(); 
		double mass=msun;
		sys.add(new Star(new double[] {4*au,0,0},new double[] {0,-5E4*mass,0},rsun,mass,5800));
	
		sys.add(new Star(new double[] {0,4*au,0},new double[] {0,0,-5E4*mass},1.7*rsun,mass,9940));
		
		sys.add(new Star(new double[]{0,0,4*au},new double[] {-5E4*mass,0,0},rsun*0.5,mass,3660));
		
	
		return sys;
	} 
	public static SolarSystem symmetric()
	{
		SolarSystem sys=new SolarSystem();
		double mass=msun;
		sys.add(new Star(new double[] {-4*au,0,0},new double[] {0,-2.5E4*mass,0},rsun,mass,10000));
		sys.add(new Star(new double[] {4*au,0,0},new double[] {0,2.5E4*mass,0},rsun,mass,10000));
		mass=8E24;
		sys.add(new Planet(new double[]{0,0,au},new double[] {0E4*mass,0E4*mass,0},new double[] {0.1,0.1,1},6000,mass,Color.green));//5
	 
		return sys;
	}
	public static SolarSystem oneSun()
	{
		SolarSystem sys=new SolarSystem();
		double mass=msun;
		sys.add(new Star(new double[] {0,0,0},new double[] {0,0,0},rsun,mass,6000));
		mass=8E24;
		sys.add(new Planet(new double[]{1.3*au,0,0},new double[] {0E4*mass,0.8E5*mass,-0E4*mass},new double[] {0.1,0.3,1},6000,mass,Color.green));//5
	 
		return sys;
	}
	
	public static SolarSystem binary()//1480
	{
		SolarSystem sys=new SolarSystem();
		double mass=msun;
		sys.add(new Star(new double[] {1.5*au,0,0},new double[] {0,mass*4.4E4,0},rsun,mass,9800));
		sys.add(new Star(new double[] {-1.5*au,0,0},new double[] {0,-mass*4.4E4,0},rsun,mass,2500));

		mass=8E24;
		sys.add(new Planet(new double[]{-2.3*au,0.15*au,0},new double[] {1E3*mass,-1.5E5*mass,-0E4*mass},new double[] {0.1,0.1,1},6000,mass,Color.green));//5
		mass=8E24;
		sys.add(new Planet(new double[]{-2.3*au,-0.15*au,0},new double[] {-1E3*mass,-1.5E5*mass,-0E4*mass},new double[] {0.1,0.5,1},6000,mass,Color.cyan));//5
	 
		return sys;
	}

	//***********************
	// print state to image
	//*********************
	private static void print(BufferedImage image,SolarSystem sys, int pl, int counter)
	{
		
		Atmosphere a=sys.planets.get(pl).atm;
		//System.out.println("atmospheresize="+a.temperature.length);
		double scale=0.0000005; 
	
		int[]	center= {image.getWidth()/2,(a.temperature[0].length+image.getHeight())/2};
				if(doble)center=new int[]{(a.temperature.length+image.getWidth())/2,(a.temperature[0].length)};//{ image.getWidth()/2,image.getHeight()/2};// 
		File file=new File(name+df.format(counter)+"."+type);
	
		//draw atmosphere 
		for(int i=0;i<a.temperature.length;i++)
			for(int j=0;j<a.temperature[0].length;j++)
				image.setRGB(i,j,a.colorcode(i,j)); 
		if(doble) {
		a=sys.planets.get(pl+1).atm;
		for(int i=0;i<a.temperature.length;i++)
			for(int j=0;j<a.temperature[0].length;j++)
				image.setRGB(i,j+a.temperature[0].length,a.colorcode(i,j)); 
		//dampen orbit;
		if(counter%4==0)	
	for(int i=a.temperature.length;i<image.getWidth();i++)
			for(int j=0;j<image.getHeight();j++) {image.setRGB(i, j, dampen(image.getRGB(i,j)));}

		}
		else
		{
			//dampen orbit;
			if(counter%3==0)
			for(int i=0;i<image.getWidth();i++)
				for(int j=a.temperature[0].length;j<image.getHeight();j++) {image.setRGB(i, j, dampen(image.getRGB(i,j)));}

		}
	
		//draw celestial bodies
		for(Star st:sys.stars)
			st.draw(image,center,scale,zBuffer,0,0);
		for(Planet planet:sys.planets)
			planet.draw(image,center,scale,zBuffer,0,0);
	
		//draw clans
		/*for(Clan cl:sys.planets.get(0).atm.clans)
		{
			cl.draw(image);
		}*/
		int[]shift= {0,0};
		//draw plants
		for(Plant cl:sys.planets.get(0).atm.plants)
		{
			cl.draw(image,shift);
		}
		//draw animals
		for(Animal cl:sys.planets.get(0).atm.animals)
		{
			cl.draw(image,shift);
		}
		if(doble) {
		shift[1]=720;
		
		for(Plant cl:sys.planets.get(1).atm.plants)
		{
			cl.draw(image,shift);
		}
		//draw animals
		for(Animal cl:sys.planets.get(1).atm.animals)
		{
			cl.draw(image,shift);
		}
		}
		//print to file
		try {
				ImageIO.write(image, type, file);
			}	catch (IOException e) {	System.out.println("IOException: Problems saving file "+name);	e.printStackTrace();}
	}
	//***********************
		// print state to image only orbit
		//*********************
		private static void print(BufferedImage image,SolarSystem sys,  int counter)
		{
			double scale=0.0000009*k; 
		
			int[]	center= {image.getWidth()/2+50,image.getHeight()/2};
			File file=new File(name+df.format(counter)+"."+type);
	
			//dampen orbit;
			if(counter%(30/steps)==0)	
		for(int i=0;i<image.getWidth();i++)
				for(int j=0;j<image.getHeight();j++) {image.setRGB(i, j, dampen(image.getRGB(i,j)));}

			
			//draw celestial bodies
			for(Star st:sys.stars)
				st.draw(image,center,scale,zBuffer,0,0);
			for(Planet planet:sys.planets)
				planet.draw(image,center,scale,zBuffer,0,0);
		
			
			
			//print to file
			try {
					ImageIO.write(image, type, file);
				}	catch (IOException e) {	System.out.println("IOException: Problems saving file "+name);	e.printStackTrace();}
		}
	//***********************
	// print state to image
	//*********************
	private static void printwrite(BufferedImage image,SolarSystem sys, int pl, int counter)
	{
			
			Atmosphere a=sys.planets.get(pl).atm;
			double scale=0.0000002; 
			int[]center={a.temperature.length/2,(a.temperature.length+a.temperature[0].length)/2};//{ image.getWidth()/2,image.getHeight()/2};// 
			File file=new File(name+df.format(counter)+"."+type);
		
			//draw atmosphere 
			for(int i=0;i<image.getWidth();i++)
				for(int j=0;j<a.temperature[0].length;j++)
					image.setRGB(i,j,a.colorcode(i,j)); 
			
			//dampen orbit;
		for(int i=0;i<image.getWidth();i++)
				for(int j=a.temperature[0].length;j<image.getHeight();j++) 
				{image.setRGB(i, j, black);}


		Writing.write(image, "Plants: "+sys.planets.get(0).atm.plants.size(),50, 1330,30,white,white);
		Writing.write(image, "Animals: "+sys.planets.get(0).atm.animals.size(),1330, 1330,30,white,white);	
			//draw celestial bodies
		/*		for(Star st:sys.stars)
				st.draw(image,center,scale,zBuffer,0,0);
			for(Planet planet:sys.planets)
				planet.draw(image,center,scale,zBuffer,0,0);
		*/
			//draw clans
			/*for(Clan cl:sys.planets.get(0).atm.clans)
			{
				cl.draw(image);
			}*/
			int[]shift= {0,0};
			//draw plants
			for(Plant cl:sys.planets.get(0).atm.plants)
			{
				cl.draw(image,shift);
			}
			//draw animals
			for(Animal cl:sys.planets.get(0).atm.animals)
			{
				cl.draw(image,shift);
			}
			//print to file
			try {
					ImageIO.write(image, type, file);
				}	catch (IOException e) {	System.out.println("IOException: Problems saving file "+name);	e.printStackTrace();}
		}
	
	//mix color towards background
	private static int dampen(int rgb)
	{ 
		double factor=0.9999999;
		Color c=new Color(rgb);
		return new Color((int)(c.getRed()*factor),(int)(c.getGreen()*factor), (int)(c.getBlue()*factor)).getRGB();
	}
	
	//update state
	public void update(double t)
	{	
		int ns=stars.size(),np=planets.size(); 
		double[][] starloc=copystarlocs(),planetlocs=copyplanetlocs();
		
		//"stars with stars..."
		for(int i=0;i<ns;i++)
		{
			double[]force=new double[dim];
			
			for(int j=0;j<ns;j++)
				if(i!=j)
				{
					add(force,relgravForce(starloc[i],starloc[j],stars.get(i).mass,stars.get(j).mass));
				}
			stars.get(i).relPull(force,t);
			
		
		}
			
		//Planets
		  for(int i=0;i<np;i++)
		  {
			  double[]force=new double[dim];
			  //	... with stars
			  for(int j=0;j<ns;j++)
			  {
				  add(force,relgravForce(planetlocs[i],starloc[j],planets.get(i).mass,stars.get(j).mass));
			  } 
			  //	... with planets
			  for(int j=0;j<0;j++) 
			  {
				  if(i!=j)
					  add(force,relgravForce(planetlocs[i],planetlocs[j],planets.get(i).mass,planets.get(j).mass));
			  }
			  planets.get(i).relPull(force,t);
			  planets.get(i).moveOn(t);
			 // log(time);
			  if(atm)
			  {
				  {planets.get(i).atm.update(stars,delt);}
				//  planets.get(i).atm.updateWind(stars,delt);
			  }
			//else 
		  }
	}
	
	public void update(double t,double[][]min)
	{	
		int ns=stars.size(),np=planets.size(); 
		double[][] starloc=copystarlocs(),planetlocs=copyplanetlocs();
		
		//"stars with stars..."
		for(int i=0;i<ns;i++)
		{
			double[]force=new double[dim];
			
			for(int j=0;j<ns;j++)
				if(i!=j)
				{
					add(force,relgravForce(starloc[i],starloc[j],stars.get(i).mass,stars.get(j).mass));
				}
			stars.get(i).relPull(force,t);
			print(stars.get(i).loc);print(stars.get(i).v);
		
		}
			
		//Planets
		  for(int i=0;i<np;i++)
		  {
			  double[]force=new double[dim];
			  //	... with stars
			  for(int j=0;j<ns;j++)
			  {
				  add(force,relgravForce(planetlocs[i],starloc[j],planets.get(i).mass,stars.get(j).mass));
				  min[i][j]=Math.min(Star.norm(Planet.subtract(planetlocs[i],starloc[j])),min[i][j]);
			  } 
			 
			  planets.get(i).relPull(force,t);
			  planets.get(i).moveOn(t);
			
			  
		  }
	}
	
	
	//write down current state for later use
	private void log(double t) 
	{
		System.out.println("Time="+t);
		System.out.println("stars:");
		for(Star st:stars)
		{
			print(st.loc);
			System.out.print(", ");
			print(st.v);
			System.out.println();
		}
	
		System.out.println("planets:");
		for(Planet pl:planets)
		{
			print(pl.loc);
			System.out.print(", ");
			print(pl.v);
			System.out.println("Live animals: "+pl.atm.animals.size());
			System.out.println("Live plants: "+pl.atm.plants.size());

		}
	}
	public static double[] relgravForce(double[] loc1,double[] loc2,double mass1,double mass2)
	{
		double[]force=new double[dim];
		double d=distance(loc1,loc2);if(d<1000)d=1000;
	//	System.out.println("d="+d);
		for(int k=0;k<dim;k++) 
		force[k]+=-g*mass1*mass2*Math.pow(d,-dim)*(loc1[k]-loc2[k]);

		return force;
	}
	double[][] copyplanetlocs() 
	{
		double[][]out=new double[planets.size()][3];
		for(int i=0;i<planets.size();i++)
			for(int j=0;j<3;j++)
				out[i][j]=planets.get(i).loc[j];
		return out;
	}
	double[][] copystarlocs() 
	{
		double[][]out=new double[stars.size()][3];
		for(int i=0;i<stars.size();i++)
			for(int j=0;j<3;j++)
				out[i][j]=stars.get(i).loc[j];
		return out;
	}
	
	
	//vector manipulation
	static void add(double[] v, double[] w) 
	{
		for(int i=0;i<dim;i++)v[i]+=w[i];
	}
	private static double distance(double[] loc1, double[] loc2) 
	{
		double out=0;
		for(int i=0;i<dim;i++)
			out+=Math.pow(loc2[i]-loc1[i], 2);
		return Math.sqrt(out);
	}
	static void print(double[]v)
	{
		System.out.print("{");
		for(int i=0;i<v.length;i++)
			System.out.print(v[i]+",");
		System.out.print("}");
	}
}
