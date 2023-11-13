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

import javax.imageio.ImageIO;

public class Atmosphere 
{
	Planet pl;
	
	static String name="hundreddays30deg", type="png";
	static DecimalFormat df=new DecimalFormat("0000");
	static double sigma=5.67/Math.pow(10, 8),	//HelmholtzBolzmann constant
			delt=3;
	double r=600;
	double c=7000;//=specific heatcapacity*density
	double[][]temperature;
	static double[]zro= {0,0,0};
	
	//************
	//Constructors
	//************
	public Atmosphere(double radius) 
	{
		r=radius;
		temperature=new double[(int) (2*Math.PI*r)][(int) (Math.PI*r)];
	}
	
	//uniform temperature t
	public Atmosphere(Planet planet,double radius,double t) 
	{
		pl=planet;
		r=radius;
		temperature=new double[(int) (2*Math.PI*r)][(int) (Math.PI*r)];
		for(int i=0;i<temperature.length;i++)
			for(int j=0;j<temperature[i].length;j++)
				temperature[i][j]=t;
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
			}
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
	private double norm(double[] v) 
	{
		return Math.sqrt(dot(v,v));
	}

	private double dot(double[] v, double[] w) 
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
