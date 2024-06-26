//****************************************
// Temperature/Star.java
// author: Non-Euclidean Dreamer
// stars with temperature, mass and radius
//****************************************

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Random;


public class Star 
{
	//Physical constants
	static double g=6*Math.pow(10,-11),
			kc=9*Math.pow(10, 9),
			c=300000*3600,
			rulescope=SolarSystem.rulescope, metric=2;
	static int dim=3,background=Color.black.getRGB();
	static int white=Color.white.getRGB(),sirius=new Color(180,199,255).getRGB(), red=Color.red.getRGB(),orange=Color.orange.getRGB(),yellow=Color.yellow.brighter().brighter().getRGB(),blue=Color.blue.brighter().getRGB();
static double count=0;
	//location attributes
	double[]loc, v;//v is actually relativistic momentum
	
	//physical attributes
	double mass, radius,temperature;
	int  color;
	
	//************
	//Constructor
	//************
	public Star(double[]l,double[]vel,double r,double m, double t) 
	{
		mass=m;
		radius=r;
		loc=l;
		v=vel;
		temperature=t;
		color=color(t);
		//color=spectrum((int)count);
	//	count+=(6*256.0)/1000;
	}
	
	//******************************************************
	// gives the star a color approx. fitting the temperature
	//*****************************************************
	public static int color(double temp)
	{
		if (temp<1700) return red;
		if (temp<3700) return orange;
		if (temp<6000) return yellow;
		if (temp<8000) return white;
		else return blue;
	}
	
	//*****************************************************
	//"relativistic pull": A time step, adjusting v and loc
	//******************************************************
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
	
	//****************************
	//draw the star onto the image
	//****************************
	public void draw(BufferedImage image, int[] center, double scale, double[][] zBuffer, int x0, int y0) 
	{
		double factor=10;//exaccerated star size
		//System.out.println("Star at i="+(center[0]+(int)(loc[0])*scale)+", j="+(center[1]+(int)(loc[1])*scale)+" with radius "+(factor*radius*scale));
		for(int i=Math.max(x0, center[0]-(int)((factor*radius-loc[0])*scale));i<center[0]+(loc[0]+factor*radius)*scale&&i<image.getWidth();i++)
		{
			double s=Math.sqrt(factor*factor*radius*radius-Math.pow((i-center[0])/scale-loc[0], 2));
			for(int j=Math.max(0, center[1]-(int)((s-loc[1])*scale));j<center[1]+(loc[1]+s)*scale&&j<image.getHeight();j++)
			//if(zBuffer[i-x0][j-y0]<loc[2])
				{
				//System.out.println("i="+i+", j="+j+", color="+color);
				//zBuffer[i-x0][j-y0]=loc[2];
				image.setRGB(i, j, color);
				
			}
		}
	}	
	
	//************
	// vector norm
	//************
	static double norm(double[] vec) 
	{
		double out=0;
		for(int i=0;i<dim;i++)
			out+=vec[i]*vec[i];
		return Math.sqrt(out);
	}
	public static int spectrum(int n)
	{
		n=n%(6*256);
		if (n<256)
			return new Color(255,n,0).getRGB();
		n-=256;
		if (n<256)
			return new Color(255-n,255,0).getRGB();
		n-=256;
		if (n<256)
			return new Color(0,255,n).getRGB();
		n-=256;
		if (n<256)
			return new Color(0,255-n,255).getRGB();
		n-=256;
		if (n<256)
			return new Color(n,0,255).getRGB();
		n-=256;
			return new Color(255,0,255-n).getRGB();
	}

}
