//**********************************************
// Temperature/Clan.java
// author: Non-Euclidean Dreamer
// Life forms on the planet spawning & evolving
//**********************************************

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.Random;

public class Clan
{
	Atmosphere atm;
	double[] comfortzone;
	double[] livzone;
	double size;
	double speed;
	double memory;
	double age=0;
	double buildingpower;
	double[][] moving;
	double[][] building;
	double[] loc;
	int color;
	static Random rand=new Random();
	public Clan(Atmosphere a, double[] l)
	{
		atm=a;
		double temp=a.temperature[(int)l[0]][(int)l[1]];
		comfortzone=new double[] {temp-1,temp+1};
		livzone=new double[] {-10,+10};
		size=1; 

		speed=0;
		memory=0;
		buildingpower=0;
		loc=l;		
		
		atm.population[(int)loc[0]][(int)loc[1]]+=1;
		color=new Color(128+rand.nextInt(128), rand.nextInt(128),rand.nextInt(128)).getRGB();
		//System.out.println("new clan: ("+loc[0]+","+loc[1]+")");
	}
	
	public void evolve(double amount)
	{
		int x=rand.nextInt(6);
		if(x==0) {comfortzone[0]-=amount/(1+age);livzone[0]+=amount/(1+age);}
		else if(x==1) {comfortzone[1]+=amount/(1+age);livzone[1]-=amount/(1+age);}
		else if(x==2) livzone[0]-=amount/(1+age);
		else if(x==3)livzone[1]+=amount/(1+age);
		else if(x==5) {size+=amount;atm.population[(int)loc[0]][(int)loc[1]]+=amount;}
		else if(x==4)speed+=amount/(1+age)*atm.r*6/atm.pl.radius;
		
		if(x<5)
		{age+=amount;
		if(rand.nextDouble()*amount>0.1)
			{
				Color col=new Color(color), nju=new Color(Math.max(0,Math.min(255, col.getRed()+rand.nextInt(5)-2)),Math.max(0,Math.min(255, col.getGreen()+rand.nextInt(5)-2)),Math.max(0,Math.min(255, col.getBlue()+rand.nextInt(5)-2)));
			}
		}
	}
	
	public Clan split()
	{
		double z=rand.nextDouble(), y=z/(1-z);
		double njusize=z*size;
		size-=njusize;
		double temp=atm.temperature[(int)loc[0]][(int)loc[1]];
		double[]live=livzone.clone(),comf=comfortzone.clone();
		double sp=speed,m=memory,bp=buildingpower;
		double aage=1.0/y, amount=aage/(1+age),yeff=y/(1+age);
		int x=rand.nextInt(7);
		if(x==0) {comf[0]+=amount; live[0]-=amount; comfortzone[0]-=yeff;livzone[0]+=yeff;}
		else if(x==1) {comfortzone[1]+=yeff;comf[1]-=amount;live[1]+=amount;livzone[1]-=yeff;}
		else if(x==2) {livzone[0]-=yeff;live[0]+=amount;}
		else if(x==3) {livzone[1]+=yeff;live[1]-=amount;}
		else if(x==4) {speed+=yeff;sp-=amount;}
		else if(x==6) {memory+=yeff;m-=amount;}
		else if(x==5) {buildingpower+=yeff;bp-=amount;}
		
		Clan out= new Clan(atm,loc.clone());
		out.buildingpower=bp;
		out.memory=m;
		out.size=njusize;
		out.speed=sp;
		if(x>4) {
		out.age=age+aage;
		age-=y;	}
		int rshift=rand.nextInt(5)-2,gshift=rand.nextInt(5)-2,bshift=rand.nextInt(5)-2;
		Color col=new Color(color);
		out.color=new Color(Math.max(0,Math.min(255, col.getRed()+rshift)),Math.max(0,Math.min(255, col.getGreen()+gshift)),Math.max(0,Math.min(255, col.getBlue()+bshift))).getRGB();
		color=	new Color(Math.max(0,Math.min(255, col.getRed()-rshift)),Math.max(0,Math.min(255, col.getGreen()-gshift)),Math.max(0,Math.min(255, col.getBlue()-bshift))).getRGB();
	
		return out;
	}
	
	public void move()
	{
		//System.out.println("speed="+speed);
		double[] gr=atm.gradient((int)loc[0],(int) loc[1]);
		double dir=((comfortzone[0]+comfortzone[1])/2-atm.temperature[(int)loc[0]][(int)loc[1]]), norm=Atmosphere.norm(gr),factor=Math.max(-speed,Math.min(dir, speed));
		if(norm==0)return;
		double j=loc[1],r=atm.r, psired=Math.acos(Math.cos(j/r)*Math.cos(factor*gr[0]/r/norm))*r+gr[1]*factor/norm,phidiff=loc[0]+r*Math.atan2(Math.sin(factor*gr[0]/r/norm),Math.sin(j/r)*Math.cos(factor*gr[0]/r/norm));
		if(psired<0) {psired=-psired; phidiff=(phidiff+atm.temperature.length/2);}
		else if(psired>=atm.temperature[0].length) {psired=2*atm.temperature[0].length-psired-.000001; phidiff=(phidiff+atm.temperature.length/2);}
		phidiff=(phidiff+atm.temperature.length)%atm.temperature.length;
		psired=Math.max(0, Math.min(atm.temperature[1].length-.01, psired));
		atm.population[(int)loc[0]][(int)loc[1]]-=size;
	//	System.out.println("move from ("+loc[0]+", "+loc[1]+") to ("+phidiff+", "+psired+")");
		loc[0]=phidiff; loc[1]=psired;
		atm.population[(int)loc[0]][(int)loc[1]]+=size;
	}

	public void draw(BufferedImage image) 
	{
		
		for(int i=(int)((loc[0]-size));i<loc[0]+size;i++)
		{
			double s=Math.sqrt(size*size-Math.pow(i-loc[0], 2));
			for(int j=-(int)((s-loc[1]));j<+(loc[1]+s);j++)
			{
				try{image.setRGB(i, j, color);}catch(ArrayIndexOutOfBoundsException e) {}
			}
		}
	}

	public double idealtemp()
	{
		return (comfortzone[1]-comfortzone[0])/2;
	}

	public void draw(BufferedImage image, int y0) 
	{
		for(int i=(int)((loc[0]-size));i<loc[0]+size;i++)
		{
			double s=Math.sqrt(size*size-Math.pow(i-loc[0], 2));
			for(int j=-(int)((s-loc[1]));j<+(loc[1]+s);j++)
			{
				try{image.setRGB(i, j+y0, color);}catch(ArrayIndexOutOfBoundsException e) {}
			}
		}
	}
}
