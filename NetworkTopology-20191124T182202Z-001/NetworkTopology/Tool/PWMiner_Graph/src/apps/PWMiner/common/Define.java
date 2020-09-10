package apps.PWMiner.common;


import java.io.File;
import java.io.IOException;

import apps.PWMiner.util.ConfigFile;

/**
 * 
 * @author ITCS 6114 QIONG CHENG
 *
 */
public class Define {
	public static String COMMAND_FILE="align.exe";
	public static String HHCOMMAND_FILE="alignHH.exe";	
	public static String HHwoCOMMAND_FILE="alignHHwo.exe";
	public static String ALIGN_RESULTS="results";
	public static String GAPLIST_FILESUFFIX=".gaplist.txt";
	public static String MISMATCHLIST_FILESUFFIX=".mismatchlist.txt";
	
	public static String DATA_DIR="data";
	public static String TRAIL_DIR="testbed";	
	public static String GRAPH_VIZ_DIR="vis";	
	public static String GRP_SUFFIX=".grp";	
	public static String VIZ_G_SUFFIX=".png";	
	public static String GRAPH_CMP_RES_DIR="results";
	public static String BATCH_MAP_DIR="mappings";
	public static String BATCH_MAP_SUFFIX=".txt";
	
	public static String COMMAND_PATH ="bin";	
	public static String WEBINF_PATH ="WEB-INF";
	//==>4 JApplet
	public static boolean isGeneralGraph=true;
	//==>
	
	//==>the information that aims to connect with url
	public static String urlConfig = "/custom.URL.config" ;
	public static String urlServiceName = "ServiceName" ;
	public static String urlHost = ".Host" ;
	public static String urlPort = ".Port" ;
	public static String urlUser = ".User" ;
	public static String urlPass = ".Password" ;
	public static String urlProtocal = "http" ;
	public static String urlFile = "/uploadgf" ;
	//==<
  
	public static File getFile(String relPathFile){
		String filename=null;
		File file = null;
		try{
			String baseDir = ConfigFile.getBaseDir(); 
			System.out.println("-----------graph file name = " + relPathFile);
			filename = baseDir + ((baseDir.lastIndexOf('\\') == baseDir.length()-1) || (baseDir.lastIndexOf('/') == baseDir.length()-1)? "" : "/" ) + DATA_DIR + "/" + relPathFile;
			System.out.println("-----------after = " + filename);
			
			file = new File(filename);
		}catch(Exception ex){
			ex.printStackTrace();
		}
		return file;
	}
	
	public static File getAbsFile(String filename){
		File file = null;
		try{			
			file = new File(filename);
		}catch(Exception ex){
			ex.printStackTrace();
		}
		return file;
	}
	
	public static boolean isGRPFile(String relPathFile){
		File file = null;
		file = getFile(relPathFile);
		
		if (relPathFile.endsWith(GRP_SUFFIX) && file != null && file.isFile()) 	return true;
		else return false;
	}
	
	public static void error(String stderr, IOException e, String where){
		System.out.println("----- Error! : " + stderr + " [" + where + "]");
		e.printStackTrace();
	}
	
	
}
