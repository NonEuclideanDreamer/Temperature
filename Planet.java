//*****************************************************************************************
// Temperature/Planet.java
// author: Non-Euclidean Dreamer
// stars with temperature, mass and radius as well as rotation information axis & daylength
//*****************************************************************************************

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.Random;

public class Planet
{
	//Physical constants
	static double g=6*Math.pow(10,-11),
			kc=9*Math.pow(10, 9),
			c=300000*3600;
	static double rulescope=SolarSystem.rulescope, metric=2;
	static int dim=3,background=Color.black.getRGB(), green=Color.green.getRGB();
		
	//location& rotation attributes
	double[]loc, v,axis,phinull,phi90;//v is the relativistic momentum, phinull, phi90 help with cartography
	double sidday,greenwich; //greenwich is the current state of rotation
	
	//physical attributes
	double mass, radius;
	int  color;
	
	public static double atmradius=300;//how fine grained we simulate the atmosphere
	
	Atmosphere atm;
	
	//***********
	//Constructor
	//***********
	public Planet(double[]l,double[]vel,double[] ax, double r,double m,Color cl) 
	{
		mass=m;
		radius=r;
		loc=l;
		v=vel;
		axis=normalize(ax);
		phinull=normalize(normcomp(ax,toplane(ax)));
		phi90=cross(axis,phinull);
		atm=new Atmosphere(this,atmradius,273);
		sidday=48;
		color=cl.getRGB();
	}

	//******************************
	// Change the Planets attributes
	//******************************
	
	//orbit on
	public void relPull(double[] force, double t) 
	{
		for(int i=0;i<dim;i++)
		{
			v[i]+=t*force[i];
		}
		double pnorm=norm(v);
		for(int i=0;i<dim;i++)
		{
			double speed=v[i]*c/Math.sqrt(Math.pow(pnorm, 2)+Math.pow(mass*c, 2));
			loc[i]+=t*speed;
		}
	}
	
	//rotate on
	public void moveOn(double t)
	{
		greenwich+=t/sidday;
	}
	
	public static void setAtmradius(double r)
	{
		atmradius=r;
	}
	
	//*********************
	// Draw planet to image
	//*********************
	public void draw(BufferedImage image, int[] center, double scale) 
	{
		double factor=1000;
		for(int i=center[0]-(int)((factor*radius-loc[0])*scale);i<center[0]+(loc[0]+factor*radius)*scale;i++)
		{
			double s=Math.sqrt(factor*factor*radius*radius-Math.pow((i-center[0])/scale-loc[0], 2));
			for(int j=center[1]-(int)((s-loc[1])*scale);j<center[1]+(loc[1]+s)*scale;j++)
			{
				image.setRGB(i, j, color);
			}
		}
	}
	
	//***************
	//vector methods
	//***************
	private double norm(double[] vec) 
	{
		double out=0;
		for(int i=0;i<dim;i++)
			out+=vec[i]*vec[i];
		return Math.sqrt(out);
	}

	public double[]normalize(double[]v)
	{
		double norm=norm(v);
		double[]out=new double[dim];
		for(int i=0;i<dim;i++)
		{
			out[i]=v[i]/norm;
		}
		return out;
	}
	//project to xy
	public double[] toplane(double[]v)
	{
		double[]out=new double[3];
		for(int i=0;i<2;i++)out[i]=v[i];
		return out;
	}
	
	public double[]cross(double[]v,double[]w)
	{
		double[]out=new double[3];
		for(int i=0;i<3;i++)
		{
			int j=(i+1)%3, k=(i+2)%3;
			out[i]=v[j]*w[k]-v[k]*w[j];
		}
		return out;
	}
	
	public double[] normcomp(double[]n,double[]v)
	{
		double parpart=dot(n,v);
		return subtract(v,times(n,parpart));
	}
	
	private double[] times(double[] v, double sc) 
	{
		double[]out=new double[dim];
		for(int i=0;i<dim;i++)
		{
			out[i]=v[i]*sc;
		}
		return out;
	}

	private double[] subtract(double[] v, double[] w)
	{
		double[]out=new double[v.length];
		for(int i=0;i<v.length;i++)
			out[i]=v[i]-w[i];
		return out;
	}
	
	private double dot(double[] v, double[] w) 
	{
		double out=0;
		for(int i=0;i<v.length;i++)out+=v[i]*w[i];
		return out;
	}
	
}