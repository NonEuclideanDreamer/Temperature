//**********************************************
// Temperature/Plant.java
// author: Non-Euclidean Dreamer
// non-moving Life forms on the planet spawning & evolving
//**********************************************


import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.Random;

public class Plant 
	{
		Atmosphere atm;
		double idealtemp;
		double[] durability;//cold/warm for plant, leaves, seed
		double[]beh;//behaviour
		double size,maxsize;
		double ejectpower;//how far off do seeds land
		double[] loc;
		double energy,minreservoir,leavefficiency;
		double seedbuffer;
		double mutability;
		int leaves;
		int color;
		boolean diploid;
		Animal parent;
		static Random rand=new Random();
		
		public Plant(Atmosphere a, double[] l,int w)//spawn a seed of haploid/diploid
		{
			parent=new Animal(a,l);
			Color c=new Color(parent.color);
			parent.color=new Color(128+c.getRed(),c.getGreen(),c.getBlue()).getRGB();
			double d=2;
			atm=a;
			mutability=rand.nextDouble();
			parent.mutability=mutability;
			idealtemp=a.temperature[(int)l[0]][(int)l[1]];
			maxsize=rand.nextDouble()*rand.nextDouble()*500;
			ejectpower=rand.nextDouble()*rand.nextDouble()*10;
			if(SolarSystem.big)ejectpower*=3;
			minreservoir=rand.nextDouble()*rand.nextDouble()*200;
			seedbuffer=rand.nextDouble()*rand.nextDouble()*10;
			durability=new double[] {rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d)};
			beh=new double[] {rand.nextDouble()*rand.nextDouble()*rand.nextDouble(),rand.nextDouble()*rand.nextDouble()*rand.nextDouble(),rand.nextDouble()*rand.nextDouble(),rand.nextDouble(),rand.nextDouble()};
			leavefficiency=rand.nextDouble()*rand.nextDouble()*10;
			energy=seedbuffer;
			size=1; 
			leaves=0;
			loc=l;		
			diploid=true;
			
			atm.vegetation[(int)loc[0]][(int)loc[1]]+=1;
			color=new Color(rand.nextInt(128),128+rand.nextInt(128), 128+rand.nextInt(128)).getRGB();
		}
		
		public Plant(Atmosphere a, double[] l)//spawn a seed
		{
			double d=2;
			atm=a;
			mutability=rand.nextDouble();
			idealtemp=a.temperature[(int)l[0]][(int)l[1]];
			maxsize=rand.nextDouble()*rand.nextDouble()*500;
			ejectpower=rand.nextDouble()*rand.nextDouble()*10;
			if(SolarSystem.big)ejectpower*=3;
			minreservoir=rand.nextDouble()*rand.nextDouble()*200;
			seedbuffer=rand.nextDouble()*rand.nextDouble()*10;
			durability=new double[] {rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d)};
			beh=new double[] {rand.nextDouble()*rand.nextDouble()*rand.nextDouble(),rand.nextDouble()*rand.nextDouble()*rand.nextDouble(),rand.nextDouble()*rand.nextDouble(),rand.nextDouble(),rand.nextDouble()};

			leavefficiency=rand.nextDouble()*rand.nextDouble()*10;
			energy=seedbuffer;
			size=1; 
			leaves=0;
			loc=l;		
			diploid=false;
			atm.vegetation[(int)loc[0]][(int)loc[1]]+=1;
			color=new Color(rand.nextInt(128),128+rand.nextInt(128), rand.nextInt(128)).getRGB();
			//System.out.println("new clan: ("+loc[0]+","+loc[1]+")");
		}
		public Plant()
		{
			
		}
		
		public Plant seed()
		{
			Plant out=new Plant();
			out.atm=atm;
			out.mutability=mutability*(1+( rand.nextDouble()-mutability)*rand.nextDouble());
			out.idealtemp=idealtemp+mutability*(rand.nextDouble()-1+rand.nextDouble())*3;
			out.maxsize=maxsize*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.ejectpower=ejectpower*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			
			out.minreservoir=minreservoir*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.seedbuffer=seedbuffer*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.diploid=false;
			out.durability=new double[6];out.beh=new double[5];
			for(int i=0;i<6;i++)
				out.durability[i]=durability[i]*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			for(int i=0;i<5;i++)
				out.beh[i]=beh[i]+( rand.nextDouble()-beh[i])*mutability*rand.nextDouble();

			out.leavefficiency=leavefficiency+( rand.nextDouble()-leavefficiency)*mutability*rand.nextDouble();
			out.energy=seedbuffer;
			energy-=seedbuffer+durability[4]+durability[5];
			out.size=0; 
			out.leaves=0;
			
			double dir=2*Math.PI*rand.nextDouble(),width=2*ejectpower/atm.r;
			double[]gr= {width*Math.cos(dir),width*Math.sin(dir)};
			double j=loc[1],r=atm.r, psired=Math.acos(Math.cos(j/r)*Math.cos(gr[0]))*r+gr[1]*r,phidiff=loc[0]+r*Math.atan2(Math.sin(gr[0]),Math.sin(j/r)*Math.cos(gr[0]));
			if(psired<0) {psired=-psired; phidiff=(phidiff+atm.temperature.length/2);}
			else if(psired>=atm.temperature[0].length) {psired=2*atm.temperature[0].length-psired-.000001; phidiff=(phidiff+atm.temperature.length/2);}
			phidiff=(phidiff+atm.temperature.length)%atm.temperature.length;
			psired=Math.max(0, Math.min(atm.temperature[1].length-.01, psired));
			out.loc=new double[] {phidiff,psired};
			
			Color col=new Color(color);
			int rshift=(int)((rand.nextDouble()-0.5)*(mutability*9)),gshift=(int)((rand.nextDouble()-0.5)*(mutability*9)),bshift=(int)((rand.nextDouble()-0.5)*(mutability*9));
			out.color=new Color(Math.max(0,Math.min(255, col.getRed()+rshift)),Math.max(0,Math.min(255, col.getGreen()+gshift)),Math.max(0,Math.min(255, col.getBlue()+bshift))).getRGB();
			System.out.println("new plant: mutability=("+out.mutability+"leafefficiency"+out.leavefficiency+",idealtemp="+out.idealtemp+", maxsize="+out.maxsize+", ejectpower="+out.ejectpower+", minreservoir"+out.minreservoir+
					" durability=("+out.durability[0]+","+out.durability[1]+","+out.durability[2]+","+out.durability[3]+","+out.durability[4]+","+out.durability[5]+")");
		
			return out;
		}
		
		public Animal baby()//for bi-generational
		{
			Animal out=new Animal();
			out.parent=this;
			out.atm=atm;
			out.mutability=mutability*(1+( rand.nextDouble()-mutability)*rand.nextDouble());
			out.leavestomach=parent.leavestomach+( rand.nextDouble()-parent.leavestomach)*mutability*rand.nextDouble();
			out.plantstomach=parent.plantstomach+( rand.nextDouble()-parent.plantstomach)*mutability*rand.nextDouble();
			out.carrionstomach=parent.carrionstomach+( rand.nextDouble()-parent.carrionstomach)*mutability*rand.nextDouble();

			out.idealtemp=idealtemp+mutability*(rand.nextDouble()-1+rand.nextDouble())*2;
			out.maxsize=maxsize*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.speed=parent.speed*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.minreservoir=minreservoir*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.seedbuffer=seedbuffer*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.beh=new double[3];
			for(int i=0;i<3;i++)
				out.beh[i]=parent.beh[i]+( rand.nextDouble()-parent.beh[i])*mutability*rand.nextDouble();
			out.diploid=true;
			out.durability=new double[6];
			for(int i=0;i<6;i++)
				out.durability[i]=durability[i]*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.energy=seedbuffer;
			energy-=seedbuffer+out.growenergy();
			out.size=1; 
			out.stomachsize=1;
			out.loc=new double[] {loc[0],loc[1]};
			atm.population[(int)loc[0]][(int)loc[1]]++;
			Color col=new Color(parent.color);
			int rshift=(int)((rand.nextDouble()-0.5)*(mutability*8)),gshift=(int)((rand.nextDouble()-0.5)*(mutability*8)),bshift=(int)((rand.nextDouble()-0.5)*(mutability*8));
			out.color=new Color(Math.max(0,Math.min(255, col.getRed()+rshift)),Math.max(0,Math.min(255, col.getGreen()+gshift)),Math.max(0,Math.min(255, col.getBlue()+bshift))).getRGB();

			System.out.println("new animal: mutability=("+out.mutability+",idealtemp="+out.idealtemp+", maxsize="+out.maxsize+", speed="+out.speed+", minreservoir"+out.minreservoir+
					", leavestomach="+out.leavestomach+", plantstomach="+out.plantstomach+", carrionstomach="+out.carrionstomach+
					" durability=("+out.durability[0]+","+out.durability[1]+","+out.durability[2]+","+out.durability[3]+","+out.durability[4]+","+out.durability[5]+")");
		
			return out;
		}
		

		public void draw(BufferedImage image,int[]shift) 
		{
			double r=1+Math.log(leaves+1);
			for(int i=(int)((loc[0]-r));i<loc[0]+r;i++)
			{
				double s=Math.sqrt(r*r-Math.pow(i-loc[0], 2));
				for(int j=-(int)((s-loc[1]));j<+(loc[1]+s);j++)
				{
					try{image.setRGB(i+shift[0], j+shift[1], color);}catch(ArrayIndexOutOfBoundsException e) {}
				}
			}
		}
		
		public double growenergy()
		{
			return(1+durability[0]+durability[1])/4;
		}
		
		public void grow()
		{
			double bc=growenergy();
			energy-=bc;
			size+=1;

		//	System.out.println("grow to"+size+" out of "+maxsize+" energy="+energy+" out of "+minreservoir);
			atm.vegetation[(int) loc[0]][(int) loc[1]]++;
		}
		
		public void growLeaf()
		{
			double bc=0.15*(durability[2]+durability[3]+1+leavefficiency);
			energy-=bc;
			leaves++;
			//System.out.println("grow a leaf");
		}


		public void draw(BufferedImage image, int y0) 
		{
			double r=1+Math.log(leaves+1);
			for(int i=(int)((loc[0]-r));i<loc[0]+r;i++)
			{
				double s=Math.sqrt(r*r-Math.pow(i-loc[0], 2));
				for(int j=-(int)((s-loc[1]));j<+(loc[1]+s);j++)
				{
					try{image.setRGB(i, j+y0, color);}catch(ArrayIndexOutOfBoundsException e) {}
				}
			}
		}
	
		
}
