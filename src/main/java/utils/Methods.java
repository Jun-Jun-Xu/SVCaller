package utils;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Arrays;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author xujun on 2022-06-23-6:11 PM
 * @project SVCaller
 */
public class Methods {

    /**
     * Get nearst min num in a int[]
     * @param value
     * @param a
     * @return int
     */
    public static int binarySearch(int value, int[] a) {
        if (value <= a[0]) { return a[0]; }
        if (value >= a[a.length - 1]) { return a[a.length - 1]; }

        int result = Arrays.binarySearch(a, value);
        if (result >= 0) { return a[result]; }

        int insertionPoint = -result - 2;
            return (a[insertionPoint] - value) < (value - a[insertionPoint - 1]) ?
                    a[insertionPoint] : a[insertionPoint - 1];

    }

    public static void sortbyColumn(int arr[][], int col) {
        // Using built-in sort function Arrays.sort
        Arrays.sort(arr, new Comparator<int[]>() {

            @Override
            // Compare values according to columns
            public int compare(final int[] entry1,
                               final int[] entry2) {

                // To sort in descending order revert
                // the '>' Operator
                if (entry1[col] > entry2[col])
                    return 1;
                else
                    return -1;
            }
        });  // End of function call sort().
    }
    /**
     * Get nearst num in a int[]
     * @param value
     * @param a
     * @return int
     */
    public static int search(int value, int[] a) {

        if(value < a[0]) {
            return a[0];
        }
        if(value > a[a.length-1]) {
            return a[a.length-1];
        }

        int lo = 0;
        int hi = a.length - 1;

        while (lo <= hi) {
            int mid = (hi + lo) / 2;

            if (value < a[mid]) {
                hi = mid - 1;
            } else if (value > a[mid]) {
                lo = mid + 1;
            } else {
                return a[mid];
            }
        }
        // lo == hi + 1
        return (a[lo] - value) < (value - a[hi]) ? a[lo] : a[hi];
    }
    /**
     * Get reverse sequence.
     * @param seq
     * @return String
     */
    public static String reverseComplement(final String seq) {
        final StringBuilder sb = new StringBuilder();
        for (int i = seq.length() - 1; i >= 0; i--) {
            final char c = seq.charAt(i);
            switch (c) {
                case 'a':
                    sb.append('t');
                    break;
                case 'A':
                    sb.append('T');
                    break;
                case 'c':
                    sb.append('g');
                    break;
                case 'C':
                    sb.append('G');
                    break;
                case 'g':
                    sb.append('c');
                    break;
                case 'G':
                    sb.append('C');
                    break;
                case 't':
                    sb.append('a');
                    break;
                case 'T':
                    sb.append('A');
                    break;
                case '-':
                    break;
                default: // n
                    sb.append(c);
                    break;
            }
        }
        return sb.toString();
    }
    /**
     * Sort a hashmap by valuee
     * @param hm
     * @return HashMao
     */
    public static HashMap<Integer, Integer> sortByValue(HashMap<Integer, Integer> hm)
    {
        // Create a list from elements of HashMap
        List<Map.Entry<Integer, Integer> > list =
                new LinkedList<Map.Entry<Integer, Integer> >(hm.entrySet());

        // Sort the list
        Collections.sort(list, new Comparator<Map.Entry<Integer, Integer> >() {
            public int compare(Map.Entry<Integer, Integer> o1,
                               Map.Entry<Integer, Integer> o2)
            {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        });

        // put data from sorted list to hashmap
        HashMap<Integer, Integer> temp = new LinkedHashMap<Integer, Integer>();
        for (Map.Entry<Integer, Integer> aa : list) {
            temp.put(aa.getKey(), aa.getValue());
        }
        return temp;
    }
    /**
     * Remove duplicates in Arraylist.
     * @param list
     * @return ArrayList
     */
    public static <T> ArrayList<T> removeDuplicates(ArrayList<T> list) {

        // Create a new LinkedHashSet
        Set<T> set = new LinkedHashSet<>();
        ArrayList list1 = new ArrayList();

        // Add the elements to set
        set.addAll(list);

        // Clear the list
//        list.clear();

        // add the elements of set
        // with no duplicates to the list
        list1.addAll(set);

        // return the list
        return list1;
    }

    /**
     * Monitor the researce usage of java.
     * @return double
     */
    private double getAverageValueByLinux() {
        try {

            long delay = 50;
            List<Double> listValues = new ArrayList<Double>();
            for (int i = 0; i < 100; i++) {
                long cput1 = getCpuT();
                Thread.sleep(delay);
                long cput2 = getCpuT();
                double cpuproc = (1000d * (cput2 - cput1)) / (double) delay;
                listValues.add(cpuproc);
            }
            listValues.remove(0);
            listValues.remove(listValues.size() - 1);
            double sum = 0.0;
            for (Double double1 : listValues) {
                sum += double1;
            }
            return sum / listValues.size();
        } catch (Exception e) {
            e.printStackTrace();
            return 0;
        }
    }
    private long getCpuT() {
        long cpuUser = 0;
        long cpuSystem = 0;
        try{
            BufferedReader reader = new BufferedReader(new FileReader("/proc/stat"));
            String line = reader.readLine();
            Pattern pattern = Pattern.compile("\\D+(\\d+)\\D+(\\d+)\\D+(\\d+)\\D+(\\d+)");
            Matcher m = pattern.matcher(line);
            if (m.find()) {
                cpuUser = Long.parseLong(m.group(1));
                cpuSystem = Long.parseLong(m.group(3));
            }
        }
        catch (Exception ex){
            ex.getStackTrace();
        }
        return cpuUser + cpuSystem;
    }
}
