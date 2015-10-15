/**
>HEADER
    Copyright (c) 2014/2015 Heewook Lee heewlee@indiana.edu

    This file is part of the GRASPER suite.

    GRASPER is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GRASPER is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GRASPER.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

import java.util.regex.*;

public class SAM implements Comparable<SAM>{
    private String samline;
    private int midpoint;
    private int effectiveLength;
    private int index;
    private SAM next;
    private short flag;
    private int start;


    public short getFlag(){
	return flag;
    }

    public int getStart(){
	return start;
    }

    public boolean isSecondary(){
        short notPrimary = 0x0100;
	if( (flag & notPrimary) == 0 )
            return false;
        return true;
    }
    
    public String getSamline(){
	return this.samline;
    }
    
    public int getMidpoint(){
	return this.midpoint;
    }
    
    public int getIndex(){
	return this.index;
    }

    public int compareTo(SAM other){
	if(this.midpoint > other.midpoint)
	    return 1;
	else if(this.midpoint == other.midpoint)
	    return 0;
	else
	    return -1;
    }
    
    public void setNext(SAM n){
	this.next = n;
    }
    
    public SAM getNext(){
	return this.next;
    }

    public int length(){
	return this.effectiveLength;
    }

    public SAM(String samline){
	this.samline = samline;
	String[] tokens = samline.split("\\t");
	this.start = Integer.parseInt(tokens[3]);
	this.effectiveLength = this.getEffectiveLengthFromCIGAR(tokens[5]);
	this.midpoint = computeMidpoint();
	
	this.flag = Short.parseShort(tokens[1]);
	this.index = -1;
	this.next = null;
    }
    
    public SAM(String samline, int i){
	this(samline);
	this.index = i;
    }
    
    private int computeMidpoint(){
	return this.start + (this.effectiveLength/2);
    }

    //we only consider M and I for the effective length of READ.                                                                                              
    private int getEffectiveLengthFromCIGAR(String cigar){
        Pattern pattern = Pattern.compile("([1-9]\\d*)([MIDNSHPX=])");
        Matcher matcher = pattern.matcher(cigar);
        int sumMI = 0;
        while(matcher.find()){
            int size = Integer.parseInt(matcher.group(1));
            char type = matcher.group(2).charAt(0);
            if(type == 'M' || type == 'I')
                sumMI += size;
        }
	return sumMI;
    }

}
