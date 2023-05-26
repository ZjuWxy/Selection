
public class Match {

	public static entry_2[] readFromMgf (String st, int len){
        
        In in = new In(st);
        entry_2[] res = new entry_2[len];
        
        
        for (int i = 1; i <= len-8; i++){
            double[] eachMZ = new double[400];
            double[] eachTIC = new double[400];
            double[] eachRT = new double[400];

            String tmp1 = in.readLine();
            String tmp2 = in.readLine();

            String s1 = in.readString();
            if(s1.substring(0,6).equals("CHARGE"))
                s1 = in.readString();
            String subS1 = s1.substring(s1.length()-10);
            if(subS1.charAt(0) < '0' || subS1.charAt(0) > '9'){
                subS1 = s1.substring(s1.length()-9);
                if(subS1.charAt(0) < '0' || subS1.charAt(0) > '9'){
                    subS1 = s1.substring(s1.length()-8);
                }
            }

            double Tmz = Double.parseDouble(subS1);

            String tmp3 = in.readLine();
            String tmp4 = in.readLine();

            int k = 0;
            while(true){
                String s_1 = in.readString();
                if(s_1.equals("END")) break;
                eachMZ[k] = Double.parseDouble(s_1);
                String s_2 = in.readString();
                eachTIC[k] = Double.parseDouble(s_2);
                String s_3 = in.readString();
                eachRT[k++] = Double.parseDouble(s_3);
            }
            String tmp5 = in.readLine();

            res[i] = new entry_2(Tmz, eachMZ, eachTIC, eachRT);
        }
        return res;
    }

    public static void main(String[] args){
    	entry_2[] res = readFromMgf("UMP-1-NEG.mgf", 1313);
    }
}