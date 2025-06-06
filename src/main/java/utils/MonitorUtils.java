package utils;

/**
 * Common sorting utils with Array
 *
 * @ author: JunXu
 * @ created: 2021-08-10 : 10:03 PM
 */

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Arrays;

public class MonitorUtils {
    /**
     * @param arg
     * get the current threads
     * @return
     */
    public static int monitor (String arg){
        int currentThreads =0;
        String command = "ps aux | grep "+arg+" | wc -l ";
        try{
            String [] cmdarry ={"/bin/bash","-c",command};
            Process p =Runtime.getRuntime().exec(cmdarry,null);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                currentThreads=(Integer.parseInt(temp.replaceAll(" ",""))-1)/2;
            }
            p.waitFor();br.close();
        }
        catch (Exception ex){
            ex.getStackTrace();
        }
        return currentThreads;
    }
    /**
     * @param user task
     * get the current threads
     * @return
     */
    public static int monitor (String user, String task){
        int currentThreads =0;
        try{
            String command = "ps aux | grep "+user+" | grep "+task+" | wc -l ";
            String [] cmdarry ={"/bin/bash/","-c",command};
            Process p =Runtime.getRuntime().exec(cmdarry,null);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                currentThreads=(Integer.parseInt(temp.replaceAll(" ",""))-1)/2;
            }
            p.waitFor();
        }
        catch (Exception ex){
            ex.getStackTrace();
        }
        return currentThreads;
    }
    /**
     * @param task sample
     * get the current process ID of a task
     * @return
     */
    public static int getProcessID (String task, String sample){
        int currentProcessID =0;
        String command = "ps aux | grep "+task+" | grep "+sample+" | awk '{print $2}' ";
        try{
            String [] cmdarry ={"/bin/bash","-c",command};
            Process p =Runtime.getRuntime().exec(cmdarry,null);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                currentProcessID=(Integer.parseInt(temp.replaceAll(" ","")));
            }
            p.waitFor();br.close();
        }
        catch (Exception ex){
            ex.getStackTrace();
        }
        return currentProcessID;
    }
    /**
     * @param processID
     * get the current process ID of a task
     * @return
     */
    public static int isFinish (int processID){
        int isFinish = 0;
        try{
            String command = "ps aux | grep "+processID+" | grep R | wc -l ";
            String [] cmdarry ={"/bin/bash/","-c",command};
            Process p =Runtime.getRuntime().exec(cmdarry,null);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                isFinish=Integer.parseInt(temp.replace(" ",""));
            }
            p.waitFor();
        }
        catch (Exception ex){
            ex.getStackTrace();
        }
        return isFinish;
    }

    /**
     * @author xujun on 2022-06-23-6:07 PM
     * @project SVCaller
     */
    public static class Methods {
        public static int binarySearch(int value, int[] a) {
            if (value <= a[0]) { return a[0]; }
            if (value >= a[a.length - 1]) { return a[a.length - 1]; }

            int result = Arrays.binarySearch(a, value);
            if (result >= 0) { return a[result]; }

            int insertionPoint = -result - 2;
            return (a[insertionPoint] - value) < (value - a[insertionPoint - 1]) ?
                    a[insertionPoint] : a[insertionPoint - 1];
        }
    }
}
