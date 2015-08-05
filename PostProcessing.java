import java.io.*;
import java.util.*;

public class PostProcessing{

    private ArrayList<Record> records;
    
    public PostProcessing(){
	this.records = new ArrayList<Record>();
    }
    
    public static void main(String [] args){
	new PostProcessing().loadRecords(args[0]);
    }
    
    public void loadRecords(String sortedRecordFile){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(sortedRecordFile));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		records.add(new Record(curline));
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	removeSingleton();
	removeMulti();
	removeWrapAround();
	//removeMulti();
	processISDra3();
	processISDra2();
	StringBuffer bf = processIS2621();
	processISUnknown();
	System.out.print(bf.toString());
	processISDra5();
	processISDra6();
	processISDra4();
	processISDra2B();

	processShortDuplicativeInsertion();

	printUnresolved();
    }
    
    public void printUnresolved(){
	System.out.println("************* UNRESOLVED ***************");
	for(int i=0; i<this.records.size();i++){
	    System.out.println(this.records.get(i)+"\n");
	}
    }

    public String pairClusters2string(ArrayList<Record> pair1, ArrayList<Record> pair2, boolean left){
	StringBuffer bf = new StringBuffer();
	for(int i=0;i<pair1.size();i++){
	    for(int j=0;j<pair2.size();j++){
		if(pair1.get(i).isInsertionPair(pair2.get(j), left)){
		    bf.append(pair1.get(i) + "\n" + pair2.get(j) + "\n\n");
		    pair1.remove(i);
		    pair2.remove(j);
		    i--;
		    j--;
		    break;
		}
	    }
	}
	if( (pair1.size() + pair2.size()) > 0){
	    bf.append("***** UNPAIRED *****\n");
	    for(int i=0;i<pair1.size();i++){
		bf.append(pair1.get(i) + "\n");
	    }
	    for(int j=0;j<pair2.size();j++){
		bf.append(pair2.get(j) + "\n");
	    }
	    bf.append("**** END OF UNPAIRED ****\n");
	}
	//bf.append("\n");
	return bf.toString();
    }

    public void pairClusters(ArrayList<Record> pair1, ArrayList<Record> pair2, boolean left){
	for(int i=0;i<pair1.size();i++){
	    for(int j=0;j<pair2.size();j++){
		if(pair1.get(i).isInsertionPair(pair2.get(j), left)){
		    System.out.println(pair1.get(i) + "\n" + pair2.get(j) + "\n");
		    pair1.remove(i);
		    pair2.remove(j);
		    i--;
		    j--;
		    break;
		}
	    }
	}
	if( (pair1.size() + pair2.size()) > 0){
	    System.out.println("***** UNPAIRED *****");
	    for(int i=0;i<pair1.size();i++){
		System.out.println(pair1.get(i));
	    }
	    for(int j=0;j<pair2.size();j++){
		System.out.println(pair2.get(j));
	    }
	    System.out.println("**** END OF UNPAIRED ****\n");
	}
	//System.out.println("\n");
    }

    public void processShortDuplicativeInsertion(){
	System.out.println("#SHORT SEGMENT DUPLICATION");

	ArrayList<Record> pair1sLeft = getCandidates(true, 146300, 146400, false, "EDGE:172:263");
	ArrayList<Record> pair2sLeft = getCandidates(true, 146300, 146400, true, "EDGE:172:263");

	ArrayList<Record> pair1sRight = getCandidates(false, 146300, 146400, false, "EDGE:172:263");
	ArrayList<Record> pair2sRight = getCandidates(false, 146300, 146400, true, "EDGE:172:263");
	
	this.pairClusters(pair1sRight, pair2sRight, false);
	this.pairClusters(pair1sLeft, pair2sLeft, true);
	System.out.println("\n");
    }


    public void processISDra3(){
	System.out.println("#DUPLICATIVE INSERTION OF ELEMENT [252945 - 254050] : EDGE 137:142 into various loci");
	getISDra3FixedIns();
	//System.out.println();

	ArrayList<Record> pair1sLeft = getCandidates(true, 253050, 253350, false, "EDGE:137:142");
	ArrayList<Record> pair2sLeft = getCandidates(true, 253650, 253950, true, "EDGE:137:142");
	
	ArrayList<Record> pair1sRight = getCandidates(false, 253050, 253350, false, "EDGE:137:142");
	ArrayList<Record> pair2sRight = getCandidates(false, 253650, 253950, true, "EDGE:137:142");
	
	this.pairClusters(pair1sRight, pair2sRight, false);
	this.pairClusters(pair1sLeft, pair2sLeft, true);
	System.out.println("\n");
    }
    
        public void getISDra3FixedIns(){
	Record pair1 = null;
	int pair1Index = -1;
	Record pair2 = null;
	int pair2Index = -1;
	
	for(int i=0;i<records.size();i++){
	    Record cur = records.get(i);
	    if(cur.regular){
		if(cur.edgeName1.equals("EDGE:137:142") && cur.edgeName2.equals("EDGE:184:240") 
		   && !cur.fwd1 && !cur.fwd2
		   && cur.doesFirstOverlapWith(253050, 253350)
		   && cur.doesSecondOverlapWith(1748200, 1748500))
		    {//pair1
			pair1 = cur;
			pair1Index = i;
			records.remove(i);
			i--;
		    }
		else if(cur.edgeName1.equals("EDGE:137:142") && cur.edgeName2.equals("EDGE:184:240") 
		   && cur.fwd1 && cur.fwd2
		   && cur.doesFirstOverlapWith(253650, 253950)
		   && cur.doesSecondOverlapWith(1747700, 1747950))
		    {//pair1
			pair2 = cur;
			pair2Index = i;
			records.remove(i);
			i--;
		    }
	    }
	}
	if(pair1Index > -1 && pair2Index >-1){
	    System.out.println("*FIXED\n" + pair1 + "\n" + pair2 + "\n");
	    //records.remove(pair1Index);
	    //records.remove(pair2Index);
	}else if(pair1Index == -1 && pair2Index > -1){
	    System.out.println("*FIXED\n**MISSING**\n" + pair2 + "\n");
	    //records.remove(pair2Index);
	}else if(pair1Index > -1 && pair2Index == -1){
	    System.out.println("*FIXED\n" + pair1 + "\n**MISSING**\n");
	    //records.remove(pair1Index);
	}else{
	    System.out.println("*FIXED\n**MISSING**\n**MISSING**\n");
	}

    }

    public void processISDra2(){
	System.out.println("#FIXED DELETION OF EDGE 144:276 (~1740bp)");
	getISDra2FixedDel();
	//System.out.println();
	
	System.out.println("#DUPLICATIVE INSERTION OF ELEMENT  [674903 - 676681] : EDGE 144:276");

	ArrayList<Record> pair1sLeft = getCandidates(true, 675050, 675350, false, "EDGE:144:276");
	ArrayList<Record> pair2sLeft = getCandidates(true, 676150, 676575, true, "EDGE:144:276");
	
	ArrayList<Record> pair1sRight = getCandidates(false, 675050, 675350, false, "EDGE:144:276");
	ArrayList<Record> pair2sRight = getCandidates(false, 676150, 676575, true, "EDGE:144:276");
	
	this.pairClusters(pair1sRight, pair2sRight, false);
	this.pairClusters(pair1sLeft, pair2sLeft, true);
	System.out.println("\n");
    }
    
        //DELETION FIXED
    public void getISDra2FixedDel(){
	for(int i=0;i<records.size();i++){
	    Record cur = records.get(i);
	    if(cur.regular){
		if(cur.doesFirstMatch(674400, 674850, true, "EDGE:145:213")
		   && cur.doesSecondMatch(676750, 677000, false, "EDGE:216:277")){
		    System.out.println(cur);
		    records.remove(i);
		    i--;
		}else if(cur.doesFirstMatch(1387250, 1387550, true, "EDGE:149:325")
			 && cur.doesSecondMatch(1389450, 1389700, false, "EDGE:168:281")){
		    System.out.println(cur);
		    records.remove(i);
		    i--;
		}else if(cur.doesFirstMatch(1611150, 1611550, true, "EDGE:192:280")
			 && cur.doesSecondMatch(1613400, 1613750, false, "EDGE:147:324")){
		    System.out.println(cur);
		    records.remove(i);
		    i--;
		}else if(cur.doesFirstMatch(1952800, 1953125, true, "EDGE:271:283")
			 && cur.doesSecondMatch(1955000, 1955300, false, "EDGE:146:313")){
		    System.out.println(cur);
		    records.remove(i);
		    i--;
		}else if(cur.doesFirstMatch(2321000, 2321350, true, "EDGE:190:282")
			 && cur.doesSecondMatch(2323250, 2323600, false, "EDGE:150:242")){
		    System.out.println(cur);
		    records.remove(i);
		    i--;
		}
	    }
	}
	System.out.println();
    }


    //111:239 IS1
    public void processISUnknown(){
	System.out.println("#DUPLICATIVE INSERTION OF ELEMENT  [846950 - 848150] : EDGE 111:239 (SINGLE COPY in the REF: ANNOTATED AS hypothetical protein)");
	getISUnknownFixedIns();
	//System.out.println();
	
	ArrayList<Record> pair1sLeft = getCandidates(true, 847025, 847475, false, "EDGE:111:239");
	ArrayList<Record> pair2sLeft = getCandidates(true, 847700, 848050, true, "EDGE:111:239");
	
	ArrayList<Record> pair1sRight = getCandidates(false, 847025, 847475, false, "EDGE:111:239");
	ArrayList<Record> pair2sRight = getCandidates(false, 847700, 848050, true, "EDGE:111:239");
	
	this.pairClusters(pair1sRight, pair2sRight, false);
	this.pairClusters(pair1sLeft, pair2sLeft, true);
	System.out.println("\n");
    }
    //112:322
    public StringBuffer processIS2621(){
	StringBuffer bf = new StringBuffer();
	bf.append("#FIXED DELETION of EDGE 112:322 (~1320bp) : EDGE 112:322");
	bf.append("\n");
	bf.append(getIS2621FixedDel());
	bf.append("\n");
	
	bf.append("#DUPLICATIVE INSERTION OF ELEMENT  [881429 - 882786] : EDGE 112:322");
	bf.append(getIS2621FixedIns());
	bf.append("\n");//System.out.println();
	
	ArrayList<Record> pair1sLeft = getCandidates(true, 881550, 881850, false, "EDGE:112:322");
	ArrayList<Record> pair2sLeft = getCandidates(true, 882400, 882650, true, "EDGE:112:322");
	
	ArrayList<Record> pair1sRight = getCandidates(false, 881550, 881850, false, "EDGE:112:322");
	ArrayList<Record> pair2sRight = getCandidates(false, 882400, 882650, true, "EDGE:112:322");
    
	bf.append(this.pairClusters2string(pair1sRight, pair2sRight, false));
	bf.append(this.pairClusters2string(pair1sLeft, pair2sLeft, true));
	
	bf.append("\n\n");

	return bf;
    }

    public void processISDra5(){
	System.out.println("#DUPLICATIVE INSERTION OF ELEMENT  [1302702 - 1303640] : EDGE 27:51");
	
	ArrayList<Record> pair1sLeft = getCandidates(true, 1302600, 1303100, false, "EDGE:27:51");
        ArrayList<Record> pair2sLeft = getCandidates(true, 1303200, 1303500, true, "EDGE:27:51");

	ArrayList<Record> pair1sRight = getCandidates(false, 1302600, 1303100, false, "EDGE:27:51");
        ArrayList<Record> pair2sRight = getCandidates(false, 1303200, 1303500, true, "EDGE:27:51");
	
	this.pairClusters(pair1sRight, pair2sRight, false);
	this.pairClusters(pair1sLeft, pair2sLeft, true);
	System.out.println("\n");
    }

    public void processISDra6(){
	System.out.println("#DUPLICATIVE INSERTION OF ELEMENT  [1538035 - 1539132] : EDGE 217:222");
	getISDra6FixedIns();
	//System.out.println();

	ArrayList<Record> pair1sLeft = getCandidates(true, 1538150, 1538450, false, "EDGE:217:222");
        ArrayList<Record> pair2sLeft = getCandidates(true, 1538650, 1539000, true, "EDGE:217:222");

	ArrayList<Record> pair1sRight = getCandidates(false, 1538150, 1538450, false, "EDGE:217:222");
        ArrayList<Record> pair2sRight = getCandidates(false, 1538650, 1539000, true, "EDGE:217:222");
	
	this.pairClusters(pair1sRight, pair2sRight, false);
	this.pairClusters(pair1sLeft, pair2sLeft, true);
	System.out.println("\n");
    }
    
    public void processISDra4(){
	System.out.println("#DUPLICATIVE INSERTION OF ELEMENT  [1455892 - 1456971] : EDGE 315:331");
	
	ArrayList<Record> pair1sLeft = getCandidates(true, 1456050, 1456350, false, "EDGE:315:331");
        ArrayList<Record> pair2sLeft = getCandidates(true, 1456550, 1456850, true, "EDGE:315:331");
	
	ArrayList<Record> pair1sRight = getCandidates(false, 1456050, 1456350, false, "EDGE:315:331");
        ArrayList<Record> pair2sRight = getCandidates(false, 1456550, 1456850, true, "EDGE:315:331");
	
	this.pairClusters(pair1sRight, pair2sRight, false);
	this.pairClusters(pair1sLeft, pair2sLeft, true);
	System.out.println("\n");
    }
    
    public void processISDra2B(){
	System.out.println("#DUPLICATIVE INSERTION OF ELEMENT [177924 - 179692] : EDGE 264:274");
	
	ArrayList<Record> pair1sLeft = getCandidates(true, 177800 , 178200, false, "EDGE:264:274");
        ArrayList<Record> pair2sLeft = getCandidates(true, 179300 , 179575, true, "EDGE:264:274");
	
	ArrayList<Record> pair1sRight = getCandidates(false, 177800 , 178200, false, "EDGE:264:274");
        ArrayList<Record> pair2sRight = getCandidates(false, 179300 , 179575, true, "EDGE:264:274");
	
	this.pairClusters(pair1sRight, pair2sRight, false);
	this.pairClusters(pair1sLeft, pair2sLeft, true);
	System.out.println("\n");
    }
    
    
    public ArrayList<Record> getCandidates(boolean left, int s, int e, boolean direction, String eName){
	ArrayList<Record> candidates = new ArrayList<Record>();
	for(int i=0; i<records.size(); i++){
	    Record cur = records.get(i);
	    if(cur.regular){
		if(left){
		    if(cur.doesFirstMatch(s, e, direction, eName)){
			candidates.add(cur);
			records.remove(i);
			i--;
		    }
		}else{
		    if(cur.doesSecondMatch(s,e,direction, eName)){
			candidates.add(cur);
			records.remove(i);
			i--;
		    }
		}
	    }
	}
	return candidates;
    }

    public void getISDra6FixedIns(){
	Record pair1 = null;
	int pair1Index = -1;
	Record pair2 = null;
	int pair2Index = -1;
	for(int i=0;i<records.size();i++){
            Record cur = records.get(i);
            if(cur.regular){
		if(cur.doesFirstMatch(697275, 697550, true, "EDGE:216:277")
		   && cur.doesSecondMatch(1538150, 1538450, false, "EDGE:217:222")){
		    pair1 = cur;
		    pair1Index = i;
		    records.remove(i);
		    i--;
		}else if(cur.doesFirstMatch(697800, 698075, false, "EDGE:216:277")
			 && cur.doesSecondMatch(1538675, 1539000, true, "EDGE:217:222")){
		    pair2 = cur;
		    pair2Index = i;
		    records.remove(i);
		    i--;
		}
	    }
	}   
	if(pair1Index > -1 && pair2Index >-1){
	    System.out.println("*FIXED\n" + pair1 + "\n" + pair2 + "\n");
	    //records.remove(pair1Index);
	    //records.remove(pair2Index);
	}else if(pair1Index == -1 && pair2Index > -1){
	    System.out.println("*FIXED\n**MISSING**\n" + pair2 + "\n");
	    //records.remove(pair2Index);
	}else if(pair1Index > -1 && pair2Index == -1){
	    System.out.println("*FIXED\n" + pair1 + "\n**MISSING**\n");
	    //records.remove(pair1Index);
	}else{
	    System.out.println("*FIXED\n**MISSING**\n**MISSING**\n");
	}
    }




    public String getIS2621FixedDel(){
	StringBuffer bf = new StringBuffer();
	Record pair1 = null;
        int pair1Index = -1;
        Record pair2 = null;
        int pair2Index = -1;

	for(int i=0;i<records.size();i++){
	    Record cur = records.get(i);
	    if(cur.regular){
		if(cur.doesFirstMatch(1637700, 1638000, true, "EDGE:147:324")
		   && cur.doesSecondMatch(1639600, 1639900, false, "EDGE:107:349")){
		    bf.append(cur + "\n\n");
		    records.remove(i);
		    i--;
		}else if(cur.doesFirstMatch(1337450, 1337750, true, "EDGE:28:109")
			 && cur.doesSecondMatch(1339350, 1339600, false, "EDGE:149:325")){
		    bf.append(cur + "\n\n");
		    records.remove(i);
		    i--;
		}else if(cur.doesFirstMatch(231925, 231975, false, "EDGE:194:330")
			 && cur.doesSecondMatch(881050, 881300, true, "EDGE:111:239")){
		    pair1 = cur;
		    records.remove(i);
		    i--;
		    pair1Index = i;
		}else if(cur.doesFirstMatch(881200, 881350, true, "EDGE:111:239")
			 && cur.doesSecondMatch(883000, 883200, false, "EDGE:185:335")){
		    pair2 = cur;
		    records.remove(i);
		    i--;
		    pair2Index = i;
		}
	    }
	}
	if(pair1Index > -1 && pair2Index >-1){
	    bf.append(pair1 + "\n" + pair2 + "\n\n");
	    //records.remove(pair1Index);
	    //records.remove(pair2Index);
	}else if(pair1Index == -1 && pair2Index > -1){
	    bf.append("**MISSING**\n" + pair2 + "\n\n");
	    //records.remove(pair2Index);
	}else if(pair1Index > -1 && pair2Index == -1){
	    bf.append(pair1 + "\n**MISSING**\n\n");
	    //records.remove(pair1Index);
	}else{
	    bf.append("**MISSING**\n**MISSING**\n\n");
	}
	
	return bf.toString();
    }

    //ArrayList<Record> pair1sLeft = getCandidates(true, 881550, 881850, false, "EDGE:112:322");
    //ArrayList<Record> pair2sLeft = getCandidates(true, 882350, 882650, true, "EDGE:112:322");
    public String getIS2621FixedIns(){
	StringBuffer bf = new StringBuffer();
	Record pair1 = null;
	int pair1Index = -1;
	Record pair2 = null;
	int pair2Index = -1;
	for(int i=0;i<records.size();i++){
            Record cur = records.get(i);
            if(cur.regular){
		if(cur.doesFirstMatch(881550, 881850, false, "EDGE:112:322")
			 && cur.doesSecondMatch(3024600, 3024900, true, "EDGE:129:290")){
		    pair1 = cur;
		    pair1Index = i;
		    records.remove(i);
		    i--;
		}else if(cur.doesFirstMatch(882400, 882650, true, "EDGE:112:322")
			 && cur.doesSecondMatch(3025050, 3025400, false, "EDGE:129:290")){
		    pair2 = cur;
		    pair2Index = i;
		    records.remove(i);
		    i--;
		}
	    }
	}
	if(pair1Index > -1 && pair2Index >-1){
	    bf.append("*FIXED\n" + pair1 + "\n" + pair2 + "\n\n");
	    //records.remove(pair1Index);
	    //records.remove(pair2Index);
	}else if(pair1Index == -1 && pair2Index > -1){
	    bf.append("*FIXED\n**MISSING**\n" + pair2 + "\n\n");
	    //records.remove(pair2Index);
	}else if(pair1Index > -1 && pair2Index == -1){
	    bf.append("*FIXED\n" + pair1 + "\n**MISSING**\n\n");
	    //records.remove(pair1Index);
	}else{
	    bf.append("*FIXED\n**MISSING**\n**MISSING**\n\n");
	}
	return bf.toString();
    }

    public void getISUnknownFixedIns(){
	Record pair1 = null;
	int pair1Index = -1;
	Record pair2 = null;
	int pair2Index = -1;
	
	for(int i=0;i<records.size();i++){
	    Record cur = records.get(i);
	    if(cur.regular){
		if(cur.edgeName1.equals("EDGE:111:239") && cur.edgeName2.equals("EDGE:72:203") 
		   && !cur.fwd1 && cur.fwd2
		   && cur.doesFirstOverlapWith(847025, 847475)
		   && cur.doesSecondOverlapWith(1820800, 1821075))
		    {//pair1
			pair1 = cur;
			pair1Index = i;
			records.remove(i);
			i--;
			
		    }
		else if(cur.edgeName1.equals("EDGE:111:239") && cur.edgeName2.equals("EDGE:72:203") 
		   && cur.fwd1 && !cur.fwd2
		   && cur.doesFirstOverlapWith(847700, 848050)
		   && cur.doesSecondOverlapWith(1821225, 1821625))
		    {//pair1
			pair2 = cur;
			pair2Index = i;
			records.remove(i);
			i--;
		    }
	    }
	}
	if(pair1Index > -1 && pair2Index >-1){
	    System.out.println("*FIXED\n" + pair1 + "\n" + pair2 + "\n");
	    //records.remove(pair1Index);
	    //records.remove(pair2Index);
	}else if(pair1Index == -1 && pair2Index > -1){
	    System.out.println("*FIXED\n**MISSING**\n" + pair2 + "\n");
	    //records.remove(pair2Index);
	}else if(pair1Index > -1 && pair2Index == -1){
	    System.out.println("*FIXED\n" + pair1 + "\n**MISSING**\n");
	    //records.remove(pair1Index);
	}else{
	    System.out.println("*FIXED\n**MISSING**\n**MISSING**\n");
	}
    }



    public void removeSingleton(){
	for(int i=0;i<records.size();i++){
	    if(records.get(i).cluID1 == 967){
		System.err.println(records.get(i));
		System.err.println(records.get(i).singleton);
	    }
	    if(records.get(i).singleton){
		System.out.println("SINGLETON:"+records.get(i).recordLine);
		records.remove(i);
		i--;
	    }
	}
    }
    
    public void removeMulti(){
	System.out.println("################## MULTI ######################");
	for(int i=0;i<records.size();i++){
	    Record cur = records.get(i);
	    if(cur.multi){
		System.out.println(cur);
		records.remove(i);
		i--;
	    }
	}
	System.out.println("############## END OF MULTI ###################");
    }

    public void removeWrapAround(){
	for(int i=0;i<records.size();i++){
	    Record cur = records.get(i);
	    if(cur.regular){
		if(cur.edgeName1.equals("EDGE:0:182") && cur.edgeName2.equals("EDGE:26:267") 
		   && !cur.fwd1 && cur.fwd2
		   && cur.doesFirstOverlapWith(100, 400)
		   && cur.doesSecondOverlapWith(2648250, 2648550))
		    {//chromosome 1
			System.out.println("WRAP-AROUND FOR chromosome 1:" + records.get(i).recordLine);
			records.remove(i);
			i--;
		    }else if(cur.edgeName1.equals("EDGE:26:267") && cur.edgeName2.equals("EDGE:129:290") 
			     && !cur.fwd1 && cur.fwd2
			     && cur.doesFirstOverlapWith(2648750, 2649050)
			     && cur.doesSecondOverlapWith(3060550, 3060900))
		    {//chromosome 2
			System.out.println("WRAP-AROUND FOR chromosome 2:" + records.get(i).recordLine);
			records.remove(i);
			i--;
		    }else if(cur.edgeName1.equals("EDGE:129:290") && cur.edgeName2.equals("EDGE:138:140") 
			     && !cur.fwd1 && cur.fwd2
			     && cur.doesFirstOverlapWith(3061100, 3061400)
			     && cur.doesSecondOverlapWith(3238050, 3238350))
		    {//MP1
			System.out.println("WRAP-AROUND FOR PLASMID MP1:" + records.get(i).recordLine);
			records.remove(i);
			i--;
		    }else if(cur.edgeName1.equals("EDGE:138:140") && cur.edgeName2.equals("EDGE:29:86") 
			     && !cur.fwd1 && cur.fwd2
			     && cur.doesFirstOverlapWith(3238600, 3238900)
			     && cur.doesSecondOverlapWith(3283500, 3283850))
		    {//CP1
			System.out.println("WRAP-AROUND FOR PLASMID CP1:" + records.get(i).recordLine);
			records.remove(i);
			i--;
		    }
	    }
	}
	System.out.println("\n");
    }

	
}

class Record{
    /*
    static boolean sortBasedOn1st = true;

    public int compareTo(Record other){
	if(sortBasedOn1st)
	    return this.start1 - other.start1;
	else
	    return this.start1 - other.start1;
	    }*/

    String recordLine;
    
    boolean singleton = false;
    boolean regular = false;
    boolean multi = false;
    
    int cluID1;
    int numEdges1;
    int numReads1;
    int mul1;
    int start1;
    int end1;
    boolean fwd1;
    String edgeName1;
    
    int cluID2;
    int numEdges2;
    int numReads2;
    int mul2;
    int start2;
    int end2;
    boolean fwd2;
    String edgeName2;

    public boolean doesFirstMatch(int s, int e, boolean dir, String eName){
	if(eName.equals(edgeName1) && (dir == fwd1) && doesFirstOverlapWith(s, e))
	    return true;
	return false;
    }

    public boolean doesSecondMatch(int s, int e, boolean dir, String eName){
	if(eName.equals(edgeName2) && (dir == fwd2) && doesSecondOverlapWith(s, e))
	    return true;
	return false;
    }
    
    public boolean doesFirstOverlapWith(int s, int e){
	if(s < start1){
	    if(e < start1)
		return false;
	    return true;
	}else{
	    if(end1 < s)
		return false;
	    return true;
	}
    }
    
    public boolean doesSecondOverlapWith(int s, int e){
	if(s < start2){
	    if(e < start2)
		return false;
	    return true;
	}else{
	    if(end2 < s)
		return false;
	    return true;
	}
    }


    public Record(String line){
	String[] tokens = line.split("\\t");
	this.recordLine = line;
	if(tokens.length == 16){
	    this.cluID1 = Integer.parseInt(tokens[0].substring(1));
	    this.numEdges1 = Integer.parseInt(tokens[1]);
	    this.numReads1 = Integer.parseInt(tokens[2]);
	    this.mul1 =Integer.parseInt(tokens[3]);
	    this.start1 = Integer.parseInt(tokens[4]);
	    this.end1 = Integer.parseInt(tokens[5]);
	    if(tokens[6].equals("--->"))
		this.fwd1 = true;
	    else if(tokens[6].equals("<---"))
		this.fwd1 = false;
	    this.edgeName1 = tokens[7];

	    this.cluID2 = Integer.parseInt(tokens[8]);
	    this.numEdges2 = Integer.parseInt(tokens[9]);
	    this.numReads2 = Integer.parseInt(tokens[10]);
	    this.mul2 =Integer.parseInt(tokens[11]);
	    this.start2 = Integer.parseInt(tokens[12]);
	    this.end2 = Integer.parseInt(tokens[13]);
	    if(tokens[14].equals("--->"))
		this.fwd2 = true;
	    else if(tokens[14].equals("<---"))
		this.fwd2 = false;
	    this.edgeName2 = tokens[15];
	    this.regular = true;
	}else if(tokens.length < 16){
	    if( (tokens[1].equals("0") && tokens[2].equals("0"))
		|| (tokens[9].equals("0") && tokens[10].equals("0")) )
		this.singleton = true;
	    else{
		System.err.println("ODD:" + line);
	    }
	}else if(tokens.length > 16){
	    this.multi = true;
	}
    }

    //if left is true --> left is the anchoring End(IS element-known)
    //--> then the right side is the insertion position so we check the rightside.
    public boolean isInsertionPair(Record other, boolean left){
	if(left){
	    //--this--><--other--
	    if(this.fwd2 && !other.fwd2){
		if( (this.start2 < other.start2)
		    && (other.end2 > this.end2)
		    && ( Math.abs(this.end2 - other.start2) <= 100 )
		    )
		    return true;
	    }
	    // --other--><--this--
	    else if(other.fwd2 && !this.fwd2){
		if( (other.start2 < this.start2)
		    && (this.end2 > other.end2)
		    && ( Math.abs(other.end2 - this.start2) <= 100 )
		    )
		    return true;
	    }
	}else{
	    //--this--><--other--
	    if(this.fwd1 && !other.fwd1){
		if( (this.start1 < other.start1)
		    && (other.end1 > this.end1)
		    && ( Math.abs(this.end1 - other.start1) <= 100 )
		    )
		    return true;
	    }
	    // --other--><--this--
	    else if(other.fwd1 && !this.fwd1){
		if( (other.start1 < this.start1)
		    && (this.end1 > other.end1)
		    && ( Math.abs(other.end1 - this.start1) <= 100 )
		    )
		    return true;
	    }
	}
	return false;
    }

    public String toString(){
	return this.recordLine;
    }
    
    
}
