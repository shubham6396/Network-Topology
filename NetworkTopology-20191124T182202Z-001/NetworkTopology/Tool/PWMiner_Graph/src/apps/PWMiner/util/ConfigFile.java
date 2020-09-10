/*
 * Created on Nov, 2007
 *
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
package apps.PWMiner.util;

import java.io.IOException;
import java.util.Properties;

import apps.PWMiner.common.*;

/**
 * @author Qiong Cheng
 *
 * To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
public class ConfigFile {

    //protected static Logger logger = Logger.getLogger(ConfigFile.class);
    
    private static String HEADLESS = "java.awt.headless";
    private static String pBaseDir = "baseDir";
    
    private static String PenalityForGap = "PenalityForGap";
    private static String PenalityForInsertion = "PenalityForInsertion";
    private static String BalanceWeight = "BalanceWeight";
    private static String EC3Cost = "EC3Cost";
    private static String EC4Cost = "EC4Cost";
    private static String ObsLogIndex = "ObsLogIndex";
    private static String PVLimit = "PVLimit";
    private static String AlignType = "AlignType";
    
    private static String ServiceName = "ServiceName";
    
    private static boolean init = false;
    static Properties props = new Properties();


    private ConfigFile() {
    }

    public static void init() {
        if (init) {
            return;
        }
        init = true;
        try {
            props.load(ConfigFile.class.getClassLoader().getResourceAsStream("/config.properties"));
        } catch (IOException e) {
            Define.error("config.properties file doesnot exist.", e, "ConfigFile.class");
        }
        if ( getProperty(HEADLESS).trim().compareToIgnoreCase("true") == 0)
        	System.setProperty(HEADLESS,"true");
        else if ( getProperty(HEADLESS).trim().compareToIgnoreCase("false") == 0)
        	System.setProperty(HEADLESS,"false");
        
        //String cdir = getBaseDir();
        //if (cdir != null ) System.setProperty("user.dir", cdir); 
        
    }

    public static String getProperty(String str) {
        init();
        return props.getProperty(str);
    }

    public static String getProperty(String str, String sdefault) {
        init();
        String value = props.getProperty(str);
        if (value != null) return value;
        else return sdefault;
    }

    
    public static String getBaseDir() {
        String basedir = getProperty(pBaseDir);
        if ( basedir.endsWith("/") || basedir.endsWith("\\")){
        	basedir = basedir.substring(0, basedir.length()-1 );
        }
        return basedir;
    }
    
    public static double getPenalityForInsertion() {
        String sPenalityForGap = getProperty(PenalityForInsertion);
        
        return new Double(sPenalityForGap).doubleValue();
    }
    
    public static double getPenalityForGap() {
        String sPenalityForGap = getProperty(PenalityForGap);
        
        return new Double(sPenalityForGap).doubleValue();
    }
    
    public static double getBalanceWeight() {
        String sBalanceWeight = getProperty(BalanceWeight);
        
        return new Double(sBalanceWeight).doubleValue();
    }
    
    public static double getEC3Cost() {
        String sEC3Cost = getProperty(EC3Cost);
        
        return new Double(sEC3Cost).doubleValue();
    }
    
    public static double getEC4Cost() {
        String sEC4Cost = getProperty(EC4Cost);
        
        return new Double(sEC4Cost).doubleValue();
    }
    
    public static double getPVLimit() {
        String sPVLimit = getProperty(PVLimit);
        
        return new Double(sPVLimit).doubleValue();
    }
    
    public static int getObsLogIndex() {
        String sObsLogIndex = getProperty(ObsLogIndex);
        
        return new Integer(sObsLogIndex).intValue();
    }
    
    public static int getAlignType() {
        String sAlignType = getProperty(AlignType);
        
        System.out.println("--------AlignType=" + sAlignType);
        return new Integer(sAlignType).intValue();
    }
    
    public static boolean isLoaded(){
    	//miner.doComparison("120", true, 0.501,2.0,10.0,1.00001,1.0);
    	//ConfigFile.getAlignType(), ConfigFile..getPenalityForGap(), ConfigFile.getBalanceWeight(), ConfigFile.getEC3Cost(), ConfigFile.getEC4Cost(), ConfigFile.getObsLogIndex(), ConfigFile.getPVLimit()
    	return init;
    }
    
    public static String getServiceName(){
    	return getProperty(ServiceName);
    }
}
