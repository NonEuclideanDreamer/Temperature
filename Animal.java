//**********************************************
// Temperature/Animal.java
// author: Non-Euclidean Dreamer
// moving Life forms on the planet spawning & evolving(feeds on plants&nutrition from diseased life forms)
//**********************************************


import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.Random;

public class Animal 
	{
		Atmosphere atm;
		double idealtemp;
		Plant parent;
		double[] durability;//cold/warm for plant, leaves, seed
		double[]beh;//behaviour
		double size,maxsize,stomachsize;
		double age=0;
		double speed;
		double[] loc;
		double energy,minreservoir,seedbuffer;
		double leavestomach,plantstomach, carrionstomach;
		double mutability;
		int color;
		boolean diploid;
		static Random rand=new Random();
		public Animal(Atmosphere a, double[] l)//spawn a seed
		{
			diploid=false;
			beh=new double[] {rand.nextDouble()*rand.nextDouble(),rand.nextDouble(),rand.nextDouble()};
			double d=2;
			atm=a;
			mutability=rand.nextDouble();
			idealtemp=a.temperature[(int)l[0]][(int)l[1]];
			maxsize=rand.nextDouble()*rand.nextDouble()*15;
			speed=rand.nextDouble()*rand.nextDouble()*50;
			if(SolarSystem.big)speed*=3;
			minreservoir=rand.nextDouble()*rand.nextDouble()*200;
			durability=new double[] {rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d),rand.nextDouble(d)};
			leavestomach=rand.nextDouble()*rand.nextDouble();
			plantstomach=rand.nextDouble()*rand.nextDouble();
			carrionstomach=rand.nextDouble()*rand.nextDouble();
			seedbuffer=rand.nextDouble()*rand.nextDouble()*10;
			energy=seedbuffer;
			size=1; 
			loc=l;		
			stomachsize=1;
			atm.population[(int)loc[0]][(int)loc[1]]+=1;
			color=new Color(rand.nextInt(128),rand.nextInt(128),128+ rand.nextInt(128)).getRGB();
			System.out.println("new animal: mutability=("+mutability+",idealtemp="+idealtemp+", maxsize="+maxsize+", speed="+speed+", minreservoir"+minreservoir+
					", leavestomach="+leavestomach+", plantstomach="+plantstomach+", carrionstomach="+carrionstomach+
					" durability=("+durability[0]+","+durability[1]+","+durability[2]+","+durability[3]+","+durability[4]+","+durability[5]+")");
		}
		public Animal()
		{
			
		}
		public Plant seedout()
		{
			Plant out=new Plant();
			out.parent=this;
			out.atm=atm;
			out.mutability=mutability*(1+( rand.nextDouble()-mutability)*rand.nextDouble());
			out.idealtemp=idealtemp+mutability*(rand.nextDouble()-1+rand.nextDouble())*2;
			out.maxsize=maxsize*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.ejectpower=parent.ejectpower*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			
			out.minreservoir=minreservoir*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.seedbuffer=seedbuffer*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.beh=new double[5];
			for(int i=0;i<5;i++)
				out.beh[i]=parent.beh[i]+( rand.nextDouble()-parent.beh[i])*mutability*rand.nextDouble();
			out.durability=new double[6];
			for(int i=0;i<6;i++)
				out.durability[i]=durability[i]*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);

			out.leavefficiency=parent.leavefficiency+( rand.nextDouble()-parent.leavefficiency)*mutability*rand.nextDouble();
			out.energy=seedbuffer;
			energy-=seedbuffer+durability[4]+durability[5];
			out.size=0; 
			out.leaves=0;
			out.diploid=diploid;
			double dir=2*Math.PI*rand.nextDouble(),width=2*parent.ejectpower/atm.r;
			double[]gr= {width*Math.cos(dir),width*Math.sin(dir)};
			double j=loc[1],r=atm.r, psired=Math.acos(Math.cos(j/r)*Math.cos(gr[0]))*r+gr[1]*r,phidiff=loc[0]+r*Math.atan2(Math.sin(gr[0]),Math.sin(j/r)*Math.cos(gr[0]));
			if(psired<0) {psired=-psired; phidiff=(phidiff+atm.temperature.length/2);}
			else if(psired>=atm.temperature[0].length) {psired=2*atm.temperature[0].length-psired-.000001; phidiff=(phidiff+atm.temperature.length/2);}
			phidiff=(phidiff+atm.temperature.length)%atm.temperature.length;
			psired=Math.max(0, Math.min(atm.temperature[1].length-.01, psired));
			out.loc=new double[] {phidiff,psired};
			
			Color col=new Color(parent.color);
			int rshift=(int)((rand.nextDouble()-0.5)*(mutability*9)),gshift=(int)((rand.nextDouble()-0.5)*(mutability*9)),bshift=(int)((rand.nextDouble()-0.5)*(mutability*9));
			out.color=new Color(Math.max(0,Math.min(255, col.getRed()+rshift)),Math.max(0,Math.min(255, col.getGreen()+gshift)),Math.max(0,Math.min(255, col.getBlue()+bshift))).getRGB();
			System.out.println("new plant: mutability=("+out.mutability+"leafefficiency"+out.leavefficiency+",idealtemp="+out.idealtemp+", maxsize="+out.maxsize+", ejectpower="+out.ejectpower+", minreservoir"+out.minreservoir+
					" durability=("+out.durability[0]+","+out.durability[1]+","+out.durability[2]+","+out.durability[3]+","+out.durability[4]+","+out.durability[5]+")");
		
			return out;
		}
		public Animal seed()
		{
			Animal out=new Animal();
			out.atm=atm;
			out.mutability=mutability*(1+( rand.nextDouble()-mutability)*rand.nextDouble());
			out.leavestomach=leavestomach+( rand.nextDouble()-leavestomach)*mutability*rand.nextDouble();
			out.plantstomach=plantstomach+( rand.nextDouble()-plantstomach)*mutability*rand.nextDouble();
			out.carrionstomach=carrionstomach+( rand.nextDouble()-carrionstomach)*mutability*rand.nextDouble();
			out.diploid=diploid;
			out.idealtemp=idealtemp+mutability*(rand.nextDouble()-1+rand.nextDouble())*2;
			out.maxsize=maxsize*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.speed=speed*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.minreservoir=minreservoir*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.seedbuffer=seedbuffer*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.beh=new double[3];
			for(int i=0;i<3;i++)
				out.beh[i]=beh[i]+( rand.nextDouble()-beh[i])*mutability*rand.nextDouble();

			out.durability=new double[6];
			for(int i=0;i<6;i++)
				out.durability[i]=durability[i]*(mutability*(rand.nextDouble()+rand.nextDouble()-1)/5+1);
			out.energy=seedbuffer;
			energy-=seedbuffer+out.growenergy();
			out.size=1; 
			out.stomachsize=1;
			out.loc=new double[] {loc[0],loc[1]};
			atm.population[(int)loc[0]][(int)loc[1]]++;
			Color col=new Color(color);
			int rshift=(int)((rand.nextDouble()-0.5)*(mutability*8)),gshift=(int)((rand.nextDouble()-0.5)*(mutability*8)),bshift=(int)((rand.nextDouble()-0.5)*(mutability*8));
			out.color=new Color(Math.max(0,Math.min(255, col.getRed()+rshift)),Math.max(0,Math.min(255, col.getGreen()+gshift)),Math.max(0,Math.min(255, col.getBlue()+bshift))).getRGB();

			System.out.println("new animal: mutability=("+out.mutability+",idealtemp="+out.idealtemp+", maxsize="+out.maxsize+", speed="+out.speed+", minreservoir"+out.minreservoir+
					", leavestomach="+out.leavestomach+", plantstomach="+out.plantstomach+", carrionstomach="+out.carrionstomach+
					" durability=("+out.durability[0]+","+out.durability[1]+","+out.durability[2]+","+out.durability[3]+","+out.durability[4]+","+out.durability[5]+")");
		
			return out;
		}
		
		

		double growenergy() {
			return (1+durability[0]+durability[1])/9;
		}
		public double livenergy(double diff, int sg)
		{
			return (size+stomachsize*(0.02*leavestomach+.01*plantstomach+0.008*carrionstomach))*(diff*diff)/Math.pow(1+durability[sg],2)/3000;
		}
		public void draw(BufferedImage image,int[]shift) 
		{
			double r=1+Math.log(size+1);
			for(int i=(int)((loc[0]-r));i<loc[0]+r;i++)
			{
				double s=Math.sqrt(r*r-Math.pow(i-loc[0], 2));
				for(int j=-(int)((s-loc[1]));j<+(loc[1]+s);j++)
				{
					try{image.setRGB(i+shift[0], j+shift[1], color);}catch(ArrayIndexOutOfBoundsException e) {}
				}
			}
		}
		

		
		public void grow()
		{
			double bc=growenergy();
			energy-=bc;
			size+=1;

		//	System.out.println("grow to"+size+" out of "+maxsize+" energy="+energy+" out of "+minreservoir);
			atm.population[(int) loc[0]][(int) loc[1]]++;
		}
	
		public void eatLeaf(Plant pl)
		{
			double counter=0;
			double bc=0.15*(1+pl.leavefficiency);
			while(pl.leaves>0&&counter<stomachsize)
			{
			pl.leaves--;
			energy+=bc*leavestomach;
			counter+=energy;
			}

			if (counter>stomachsize)stomachsize*=1.01;
		}
		
		public void eatPlant(Plant pl)
		{
		
			double bc=Math.min( pl.size,stomachsize);	
			pl.size-=bc;
			if (bc==stomachsize)stomachsize*=1.01;
			energy+=bc*plantstomach;
		}
		public void eatBattery(Plant pl)
		{
			double bc=Math.min(pl.energy, stomachsize);
			if (bc==stomachsize)stomachsize*=1.01;
			energy+=bc;
			pl.energy-=bc;
		}
		public void eatNutrition()
		{
			int x=(int)loc[0], y=(int)loc[1];
			double bc=Math.min(atm.nutrition[x][y], stomachsize);
			energy+=bc*carrionstomach;
			if (bc==stomachsize)stomachsize*=1.01;
			atm.nutrition[x][y]-=bc;
		}
		public void moveTemp()
		{
			//System.out.println("speed="+speed);
			int g=0;
			
			double[] gr=atm.gradient((int)loc[0],(int) loc[1]);
			double dir=(idealtemp-atm.temperature[(int)loc[0]][(int)loc[1]]), norm=Atmosphere.norm(gr),factor=Math.max(-speed,Math.min(dir, speed));
			if (dir>0)g=1;
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
			energy-=Math.abs(factor/norm*dir/durability[2+g]/30);
		}
		public void movetoFood()
		{
			int g=0;
			
			int x=(int)loc[0],y=(int) loc[1];
			double[]gr=add(gradient(atm.vegetation,x,y),gradient(atm.nutrition,x,y));
			double dir=(idealtemp-atm.temperature[(int)loc[0]][(int)loc[1]]), norm=Atmosphere.norm(gr),factor=Math.max(-speed,Math.min(dir, speed));
			if (dir>0)g=1;
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

			energy-=Math.abs(factor/norm*dir/durability[2+g]/30);
		}
		public void draw(BufferedImage image, int y0) 
		{
			double r=1+Math.log(size+1);
			for(int i=(int)((loc[0]-r));i<loc[0]+r;i++)
			{
				double s=Math.sqrt(r*r-Math.pow(i-loc[0], 2));
				for(int j=-(int)((s-loc[1]));j<+(loc[1]+s);j++)
				{
					try{image.setRGB(i, j+y0, color);}catch(ArrayIndexOutOfBoundsException e) {}
				}
			}
		}
		public double[] gradient(double[][]field,int i,int j)
		{
			int jp=Math.min(j+1,field[0].length-1),jm=Math.max(j-1,0),im=i-1,ip=i+1;
			if(im<0)im=field.length-1;
			if(ip==field.length)ip=0;
			
			return new double[] {(field[ip][j]-field[im][j])/Math.sin((j+0.1)/atm.r),field[i][jp]-field[i][jm]};
		}
		public double[] add(double[]v,double[]w)
		{
			double[]out=new double[v.length];
			for(int i=0;i<v.length;i++)
			{
			out[i]=v[i]+w[i]	;
			}
			return out;
		}
		public void eatSeed(Plant pl) 
		{
			energy+=pl.energy*plantstomach*plantstomach;
			pl.energy=0;
		}
}
