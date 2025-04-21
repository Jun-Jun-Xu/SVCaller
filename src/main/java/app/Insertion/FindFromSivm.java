package app.Insertion;

import java.lang.reflect.Parameter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @author xujun on 2022-11-26-5:17 PM
 * @project SVCaller
 */
public class FindFromSivm {
    HashSet queHS = null;
    String inputFile = null;
    String paraFile = null;
    String outputFile = null;
    HashSet[] nameHS = null; HashMap<String, Integer>[] nameFreqHM = null;
    HashMap<String, ArrayList>[] nameInfor = null;
    ArrayList[] result = null; ArrayList[] resultS = null;

    public FindFromSivm(String[] parameters) {
        long startTime = System.currentTimeMillis();
        queHS = new HashSet();
        queHS.add("M"); queHS.add("I"); queHS.add("S"); queHS.add("="); queHS.add("X");
        inputFile = parameters[0];
        paraFile = parameters[1];
        outputFile = parameters[2];
        this.running();
        long endTime   = System.currentTimeMillis();
        NumberFormat formatter = new DecimalFormat("#0.00000");
        System.out.println(formatter.format((endTime - startTime) / 1000d));
    }
    public void running(){

    }
}
