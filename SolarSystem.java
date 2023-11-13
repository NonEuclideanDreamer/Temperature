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
	
	static String name="test",type="png";
	static DecimalFormat df=new DecimalFormat("0000");
	
	//physical constants in km,kg,h units
	static double g=8.6E-13;
	static int dim=3; 
	static double rulescope=dim-1;
	static double[]zro= {0,0,0};
	static double metric=2, delt=4;
	static int steps=5,black=Color.black.getRGB();//steps:How many iteration steps "+delt" between images	
	static double rsun=690000,msun=2E30,au=1.5E8;

	static Random rand=new Random();
	
	//*************
	// Constructor
	//*************
	public SolarSystem() 
	{
		planets=new ArrayList<Planet>(); 
		stars=new ArrayList<Star>();
	}
	private void add(Planet planet) 
	{
		planets.add(planet);	
	}

	private void add(Star star) 
	{
		stars.add(star);
	}
	
	//*****************
	// Main Method
	//************
	public static void main(String[]args) 
	{
		SolarSystem sys=wild();

		int counter=0;
		Atmosphere a=sys.planets.get(0).atm;
		BufferedImage image=new BufferedImage(a.temperature.length,a.temperature.length,BufferedImage.TYPE_4BYTE_ABGR);
		while(true)  
		{ 
			System.out.println(counter); 
			for(int i=0;i<steps;i++)	sys.update(delt);
			print(image,sys,0,counter);
			counter++; 
		}
	}
	
	//************************
	// Specific Solar Systems
	//************************
	public static SolarSystem stablish()
	{
		SolarSystem sys=new SolarSystem();
		double mass=2*msun;
		sys.add(new Star(new double[] {-4*au,0,0},new double[] {0,-4.2E4*mass,0},rsun,mass,5500));
		sys.add(new Star(new double[] {4*au,0,0},new double[] {0,4.2E4*mass,0},rsun,mass,5500));
		mass=.5*msun;
		sys.add(new Star(new double[]{0,0,2*au},new double[] {0,0,0},rsun*Math.pow(0.25,1/3.0),mass,8000));
		
		mass=8E24;
		sys.add(new Planet(new double[]{3*au,0,0},new double[] {0,-1E5*mass,0},new double[] {0.1,0.1,1},6000,mass,Color.green));

		return sys;
	}
	
	public static SolarSystem wild()
	{
		SolarSystem sys=new SolarSystem();
		double mass=msun;
		sys.add(new Star(new double[] {-4*au,0,0},new double[] {0,-4.2E4*mass,0},rsun,mass,5800));
		mass=2*msun;
		sys.add(new Star(new double[] {4*au,0,0},new double[] {0,2.6E4*mass,0},1.7*rsun,mass,9940));
		
		mass=.5*msun;
		sys.add(new Star(new double[]{0,0,2*au},new double[] {0,-2E4*mass,0},rsun*0.5,mass,3660));
		
		mass=8E24;
		sys.add(new Planet(new double[]{-3*au,0,0},new double[] {0,7.5E4*mass,0},new double[] {0.1,0.1,1},6000,mass,Color.green));

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

	//***********************
	// print state to image
	//*********************
	private static void print(BufferedImage image,SolarSystem sys, int pl, int counter)
	{
		Atmosphere a=sys.planets.get(pl).atm;
		double scale=0.0000005;
		int[]center= {a.temperature.length/2,(a.temperature.length+a.temperature[0].length)/2};
		File file=new File(name+df.format(counter)+"."+type);
		
		//draw atmosphere
		for(int i=0;i<image.getWidth();i++)
			for(int j=0;j<a.temperature[0].length;j++)
				image.setRGB(i,j,a.colorcode(i,j));
		
		//dampen orbit;
		for(int i=0;i<image.getWidth();i++)
			for(int j=a.temperature[0].length;j<image.getHeight();j++)image.setRGB(i, j, dampen(image.getRGB(i,j)));

		for(Star st:sys.stars)
			st.draw(image,center,scale);
		for(Planet planet:sys.planets)
			planet.draw(image,center,scale);

		//print to file
		try {
				ImageIO.write(image, type, file);
			}	catch (IOException e) {	System.out.println("IOException: Problems saving file "+name);	e.printStackTrace();}
	}
	
	//mix color towards background
	private static int dampen(int rgb)
	{
		double factor=0.99;
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
			  for(int j=0;j<np;j++) 
			  {
				  if(i!=j)
					  add(force,relgravForce(planetlocs[i],planetlocs[j],planets.get(i).mass,planets.get(j).mass));
			  }
			  planets.get(i).relPull(force,t);
			  planets.get(i).moveOn(t);
			  planets.get(i).atm.update(stars,delt);
		  }
	}
	public static double[] relgravForce(double[] loc1,double[] loc2,double mass1,double mass2)
	{
		double[]force=new double[dim];
		double d=distance(loc1,loc2);
	//	System.out.println("d="+d);
		for(int k=0;k<dim;k++) 
		force[k]+=-g*mass1*mass2*Math.pow(d,-dim)*(loc1[k]-loc2[k]);

		return force;
	}
	private double[][] copyplanetlocs() 
	{
		double[][]out=new double[planets.size()][3];
		for(int i=0;i<planets.size();i++)
			for(int j=0;j<3;j++)
				out[i][j]=planets.get(i).loc[j];
		return out;
	}
	private double[][] copystarlocs() 
	{
		double[][]out=new double[stars.size()][3];
		for(int i=0;i<stars.size();i++)
			for(int j=0;j<3;j++)
				out[i][j]=stars.get(i).loc[j];
		return out;
	}
	
	
	//vector manipulation
	private void add(double[] v, double[] w) 
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

}
