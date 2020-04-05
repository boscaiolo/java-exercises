package org.boscaiolo.workout;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Main {

    public static int binaryGap(int n) {
        String binaryForm = Integer.toBinaryString(n);
        String[] zeros = binaryForm.split("1+");
        String[] copy = binaryForm.endsWith("0") ? Arrays.copyOf(zeros, zeros.length - 1) : zeros;
        int result = 0;
        for (String string : copy) {
            if (string.length() > result) {
                result = string.length();
            }
        }
        return result;
    }

    public static int[] rotateRight(int[] A, int k) {
        if (A == null || A.length == 0)
            return A;
        int[] a = A.clone();
        for (int j = 0; j < k; j++) {
            int last = a[a.length - 1];
            for (int i = a.length - 1; i > 0; i--) {
                a[i] = a[i - 1];
            }
            a[0] = last;
        }
        return a;
    }

    public static int orphan(int[] a) {
        int x = 0;

        for (int i : a) {
            x = x ^ i;
        }

        return x;
    }

    public static int mindiff(int[] A) {
		int total = 0;

		for (int i = 0; i < A.length; i++)
			total += A[i];

		int diff = Integer.MAX_VALUE;
		int prev = 0;
		int next = total;

		for (int p = 1; p < A.length; p++) {
			prev += A[p - 1];
			next = total - prev;
			diff = Math.min(diff, Math.abs(prev - next));
		}

		return diff;
    }

    public static int missing(int[] A){
        long sum = 0;
        
        for (int i : A) {
            sum+=i;
        }
        long n = A.length + 1;

        long f = (n*(n+1))/2;
        
        return (int)(f - sum);
    }

    public static int frogImp(int X, int Y, int D){
		if ((Y - X) % D == 0) {
			return (Y - X) / D;
		}
		
		return (Y - X) / D + 1;
    }

    public static int frogRiver(int[] A, int X){
        Set<Integer> path = IntStream.rangeClosed(1, X).parallel().boxed().collect(Collectors.toSet());

        for(int i = 0; i<A.length; i++){
            path.remove(A[i]);
            if(path.isEmpty()){
                return i;
            }
        }
        
        return -1;
    }

    public static int[] counters(int[] A, int N){
		int max = 0;
		int lastMax = 0;
		int[] counters = new int[N];

		for (int k = 0; k < A.length; k++) {
			if (A[k] >= 1 && A[k] <= N) {
				// this is instead of the nested loop in solution2
				counters[A[k] - 1] = Math.max(counters[A[k] - 1], lastMax);
				counters[A[k] - 1] += 1;
				max = Math.max(max, counters[A[k] - 1]);
			} else if (A[k] == N + 1) {
				lastMax = max;
			}
		}
		
		for(int i = 0 ; i < counters.length; i++) {
			counters[i] = Math.max(counters[i], lastMax);
		}

		return counters;
    }

    public static int smallestPosMissing(int[] A){
        Set<Integer> set = new HashSet<>();
        
        for (int i : A) {
            if(i>0){
                set.add(i);
            }
        }

        for(int i=1; i<=set.size(); i++){
            if(!set.contains(i)){
                return i;        
            }
        }
        
        return set.size() + 1;
    }

    public static int countDiv(int A, int B, int K){
        return (B/K) - (A/K) + (A % K == 0 ? 1 : 0);
    }

    public static int[] genRange(String S, int[] P, int[] Q){
        //used jagged array to hold the prefix sums of each A, C and G genoms
        //we don't need to get prefix sums of T, you will see why.
        int[][] genoms = new int[3][S.length()+1];
        //if the char is found in the index i, then we set it to be 1 else they are 0
        // 3 short values are needed for this reason
        short a, c, g;
        for (int i=0; i<S.length(); i++) {
            a = 0; c = 0; g = 0;
            if ('A' == (S.charAt(i))) {
                a=1;
            }
            if ('C' == (S.charAt(i))) {
                c=1;
            }
            if ('G' == (S.charAt(i))) {
                g=1;
            }
            //here we calculate prefix sums. To learn what's prefix sums look at here https://codility.com/media/train/3-PrefixSums.pdf
            genoms[0][i+1] = genoms[0][i] + a;
            genoms[1][i+1] = genoms[1][i] + c;
            genoms[2][i+1] = genoms[2][i] + g;
        }

        int[] result = new int[P.length];
        //here we go through the provided P[] and Q[] arrays as intervals
        for (int i=0; i<P.length; i++) {
            int fromIndex = P[i]+1;
            int toIndex = Q[i]+1;
			//if the substring contains a, then genoms[0][toIndex] - genoms[0][fromIndex-1] > 0
            if (genoms[0][toIndex] - genoms[0][fromIndex-1] > 0) {
                result[i] = 1;
            } else if (genoms[1][toIndex] - genoms[1][fromIndex-1] > 0) {
                result[i] = 2;
            } else if (genoms[2][toIndex] - genoms[2][fromIndex-1] > 0) {
                result[i] = 3;
            } else {
                result[i] = 4;
            }
        }

        return result;
    }

    public static int minAvgTwoSlice(int[] A) {
        final int N = A.length;

        int minIndex = 0;
        double minAvg = Double.MAX_VALUE;
  
        for ( int i = 0; i < N - 1; i++ ) {
          double average = ( A[ i ] + A[ i + 1 ] ) / 2.0;
  
          if ( i < N - 2 ) {
            double threeSliceAvg = ( A[ i ] + A[ i + 1 ] + A[ i + 2 ] ) / 3.0;
            average = Math.min( average, threeSliceAvg );
          }
  
          if ( average < minAvg ) {
            minAvg = average;
            minIndex = i;
          }
        }
  
        return minIndex;
    }

    public int passingCars(int[] A) {
        
        
        int countZero= 0;
        int pair = 0;
        int counter = 0;
        for (int i = A.length-1; i>=0;i--){
            counter++;
            if (A[i]==0){
                countZero++;
                pair += counter-countZero;
                if (pair > 1000000000){
                    return -1;    
                }
            }    
        }
        
        return pair;
        
    }

    public static int dist(int[] A){
        Set<Integer> s = new HashSet<>();
        for(int i : A){
            s.add(i);
        }
        return s.size();
    }

    public static int maxthree(int[] A){
        Arrays.sort(A);

        int max_1 = A[A.length-1] * A[A.length-2] * A[A.length-3];
        
        int max_2 = A[0] * A[1] * A[A.length-1];

        int max = Math.max(max_1, max_2);
        
        return max;
    }

    public static int intersectingDiscs(int[] A){
        int n = A.length;
        int[] sum = new int[n];
        
        for (int i = 0; i < n; i++) {
            int right;
            //if i+A[i]<= n-1, that's it, extract this i+A[i], let sum[i+A[i]]++, means there is one disk that i+A[i]
            if (n - i - 1 >= A[i]){
                right = i + A[i];
            } else {
                right = n - 1;    
            }
            
            sum[right]++;
        }
        
        for (int i = 1; i < n; i++) {
            sum[i] += sum[i - 1];  //sum[i] means that there are sum[i] number of values that <= i;
        }
        
        long ans = (long) n * (n - 1) / 2;
        
        for (int i = 0; i < n; i++) {
            int left;
            
            if (A[i] > i) {
                left = 0;
            } else {
                left = i - A[i];// Find the positive i-A[i].     
            }
            
            if (left > 0){
                ans -= sum[left - 1];//Find the number that is smaller than 1-A[i], sum[n-1] will never be used as we only need sum[n-1-1] at most.  
            } 
        }
        
        if (ans > 10000000) {
            return -1;    
        }
        
        return (int) ans;
    }

    public int solution(int[] A) {
         Arrays.sort(A);  
         if (A.length<3) {
            return 0;    
         }
         for (int i=0;i<A.length-2;i++) {
            if (A[i]>A[i+2]-A[i+1]){  
                return 1;
            }         
         }
         
         return 0;
    }

    public static int brackets(String S) {

        // main idea: use "Stack" (push and pop)
        
        //special case
        if(S.length() == 0)
            return 1;
        
        // new Stack<Character>()
        Stack<Character> stack = new Stack<>();
        
        // scan the string (just one pass)
        for(int i=0; i< S.length(); i++){    
            // note: push "its pair"
            if( S.charAt(i) == '(' ){
                stack.push(')');
            }
            else if( S.charAt(i) == '[' ){
                stack.push(']');
            }
            else if( S.charAt(i) == '{' ){
                stack.push('}');
            }
            // pop and check
            else if( S.charAt(i) == ')' || S.charAt(i) == ']' || S.charAt(i) == '}'){
                // important: check if the stack is empty or not (be careful)
                if(stack.isEmpty() == true){
                    return 0;
                }
                else{
                    char temp = stack.pop(); // check if the stack is empty before pop!!!
                    if(temp != S.charAt(i)){ // not a pair
                        return 0;
                    }
                }
            }
        }
        // note: check if the stack is empty or not (be careful)
        if( !stack.isEmpty() ){
            return 0;
        }
        else{
            return 1;
        }
    }

    public static int fish(int[] A , int[] B){
        // special case: no fish
        if(A.length == 0)
            return 0;

        // main idea: use "stack" to store the fishes with B[i]==1 
        // that is, "push" the downstream fishes into "stack"
        // note: "push" the Size of the downstream fish
        Stack<Integer> st = new Stack<>(); 
        int numAlive = A.length;
        
        for(int i=0; i<A.length; i++){
            
            // case 1; for the fish going to downstrem
            // push the fish to "stack", so we can keep them from the "last" one
            if(B[i] ==1){
                st.push(A[i]); // push the size of the downstream fish
            }
            // case 2: for the fish going upstream
            // check if there is any fish going to downstream
            if(B[i] ==0){
                while( !st.isEmpty() ){ 
                    // if the downstream fish is bigger (eat the upstream fish)
                    if( st.peek() > A[i] ){
                        numAlive--;
                        break;      // the upstream fish is eaten (ending)   
                    }
                    // if the downstream fish is smaller (eat the downstream fish)
                    else if(st.peek() < A[i]){
                        numAlive--;
                        st.pop();   // the downstream fish is eaten (not ending)
                    }
                }
            }  
        }    
        
        return numAlive;
    }

    public static int nested(String S){
        Stack<Character> stack = new Stack<>();

        for (char c : S.toCharArray()) {
            if(c == '('){
                stack.push(c);
            }
            if ( c == ')'){
                if(stack.isEmpty()){
                    return 0;
                } else {
                    stack.pop();
                }
            }
        }

        if (stack.isEmpty()){
            return 1;
        }
        return 0;
    }

    public static int wall(int[] A){
    
        Stack<Integer> st = new Stack<>();
        int numBlock =0;
    
        for(int i=0; i< A.length; i++){
        
            while(!st.isEmpty() && st.peek() > A[i] ){
                st.pop();
            }
            if(st.isEmpty() || st.peek() < A[i] ){
                numBlock++;    
                st.push(A[i]); 
            }
        }
        
        return numBlock;
    }

    public static int domin(int[] A){
        Map<Integer, Integer> map = new HashMap<>();
        if(A.length==1){
            return 0;
        }
        for (int i = 0; i < A.length; i++) {
                Integer k = map.get(A[i]);
                if(k==null){
                    map.put(A[i],1);
                }else{
                    k++;
                    if(k > A.length/2){
                        return i;
                    }
                    map.put(A[i], k);
                }
        }

        return -1;
    }

    public static int leaderCount(int[] A){
        if(A.length==1) {
            return 0;    
        }
        
        int value = A[0];
        int size=0;
        for(int i=0;i<A.length;i++) {
            if(size==0) {
                size++;    
                value = A[i];
            }else {
                if(A[i]==value) {
                    size++;    
                }else {
                    size--;
                }                
            }   
        }
        int candidate = -1;
        int count = 0;     
        if(size>0) {
           candidate = value;     
        }
        
        for(int i=0;i<A.length;i++) {
            if(A[i]==candidate) {
                count++;
            }    
        }

        if(count<=A.length/2) {  
            return 0;
        }
        
        int leader = candidate;
        int equiCount = 0;
        int leaderCount = 0;
        for(int i=0;i<A.length;i++) {
            if (A[i] == leader) {
                leaderCount++;    
            }
            
            if(leaderCount>(i+1)/2  && (count-leaderCount)>(A.length-i-1)/2) {
                equiCount++;    
            }
        }
        
        return equiCount;
    }
    

    public static int maxDSsum(int[] A){
        int N = A.length;
        int[] K1 = new int[N];
        int[] K2 = new int[N];

        for(int i = 1; i < N-1; i++){
          K1[i] = Math.max(K1[i-1] + A[i], 0);
        }
        for(int i = N-2; i > 0; i--){
          K2[i] = Math.max(K2[i+1]+A[i], 0);
        }

        int max = 0;

        for(int i = 1; i < N-1; i++){
          max = Math.max(max, K1[i-1]+K2[i+1]);
        }

        return max;
    }

    public static int proffit(int[] A){

		int maxProfit = 0;
		int minPrice = Integer.MAX_VALUE;
		for (int price : A) {
			minPrice = Math.min(minPrice, price);
			maxProfit = Math.max(maxProfit, price - minPrice);
		}
		return maxProfit;
    }

    public static int maxSsum(int[] A){
        // initial setting A[0]
        int maxEndingPrevious = A[0];
        int maxEndingHere = A[0];
        int maxSoFar = A[0];
 
        // note: for i=0, it will return A[0] (also for "one element" cases)  
           
        for(int i = 1; i < A.length; i++){
            maxEndingHere = Math.max(A[i], maxEndingPrevious + A[i]); // <--- key point~!!
            maxEndingPrevious = maxEndingHere;
            maxSoFar = Math.max(maxSoFar, maxEndingHere); // update the max (be careful)
        }
        
        return maxSoFar; // can be used for "all negative" cases 
    }


    public static int divisorCount(int N){
        int count = 0;

        int lim = (int) Math.sqrt(N);

        for (int i = 1; i <= lim; i++) {
            if(N%i==0){   
                count+=2;
            }        
        }        
        if(lim * lim == N){
            count--;
        }
        return count;
    }

    public static int flags(int[] A) {
        List<Integer> list = new ArrayList<>();
        for(int i = 1; i < A.length -1; i++){
        if(A[i - 1] < A[i] && A[i] > A[i+1]){
            list.add(i);
            
        }
        }
        if(list.size() == 1 || list.size() == 0){
            return list.size();
        }
        int sf = 1;  
        int ef = list.size();  
        int result = 1;  
        while (sf <= ef) {  
            int flag = (sf + ef) / 2;  
            boolean suc = false;  
            int used = 0;  
            int mark = list.get(0);  
            for (int i = 0; i < list.size(); i++) {  
                if (list.get(i) >= mark) {  
                    used++;  
                    mark = list.get(i) + flag;  
					if (used == flag) {                       
						suc = true;  
						break;  
					}  
                }  
            }  
            if (suc) {  
                result = flag;  
                sf = flag + 1;  
            }else {  
                ef = flag - 1;  
            }  
        }  
        return result; 
    }

    public static int minRect(int N) {
        int sideA = 1;
        int sideB = N;
        
        int p = (sideA + sideB)*2;

        while(sideA < sideB){
            sideA++;
            sideB = N /sideA;
            if(sideB * sideA==N && p>(sideA + sideB)*2){
                p=(sideA + sideB)*2;
            }
        }

        return p;
    }

    public static int blocks(int[] A) {
        List<Integer> list = new ArrayList<>();
        for(int i = 1; i < A.length -1; i++){
            if(A[i - 1] < A[i] && A[i] > A[i+1]){
                list.add(i);
            }
        }
        
        if(list.size() == 1 || list.size() == 0){
            return list.size();
        }
        
        int blocks = list.size();
        
        while(blocks > 1){
            if(A.length%blocks==0){
                int size = A.length/blocks;
                
                int j = 0;

                for(int peaksIndex : list){
                    if( peaksIndex/size == j){ 
                        j++; 
                    }
                }

                if(j == blocks){
                    return blocks;
                }
            }
            blocks--; 
        }

        return blocks;
    }

    public static int[] countNonDiv(int[] A){
        // main idea:
        // using "HashMap" to count 
        
        // map1(key, value)
        HashMap<Integer, Integer> map1 = new HashMap<>();
        // key: the elements, value, count of elements
        for(int i=0; i< A.length; i++){
            if(map1.containsKey(A[i]) == false){
                map1.put(A[i], 1); // add new element
            }
            else{
                 map1.put(A[i], map1.get(A[i])+1 ); // count++
            }
        }
        
        // map2(key, value)
        HashMap<Integer, Integer> map2 = new HashMap<>();
        // key: the elements, value, count of "number of non-divisors" of elements
        for( int n : map1.keySet() ){            
            int numDivisors =0;
            // find divisors from 1 to sqrt(n)
            int sqrtN = (int)Math.sqrt(n);  
            for(int i=1; i<=sqrtN; i++ ){
                if( n % i == 0){ // means: i could be a divisor
                    int anotherDivisor = n/i; 
    
                    if(map1.containsKey(i) == true ){
                        numDivisors = numDivisors + map1.get(i);
                    }
                    if(anotherDivisor != i){ // avoid double count (be careful)
                        if(map1.containsKey(anotherDivisor) == true){
                            numDivisors = numDivisors + map1.get(anotherDivisor);
                        }
                    }
                }
            }
            
            int numNonDivisors = A.length - numDivisors;
            map2.put(n, numNonDivisors); 
        }
        
        // results: number of non-divisors
        int[] results = new int[A.length];
        for (int i = 0; i < A.length; i++) {
            results[i] = map2.get(A[i]);
        }

        return results;
    }

    public static int[] semiPrimes(int N, int[] P, int[] Q){
        int[] factors = new int[N+1];
        for(int i = 2; i*i<=N; i++){
            if(factors[i] == 0){
                for(int k = i * i; k <= N; k+=i){
                        factors[k]=i;
                }
            }
        }

        int[] semiPrimes = new int[N+1];
        for(int i = 4; i<factors.length; i++){
            if(factors[i] != 0){
                int ff = factors[i];
                if(factors[i/ff]==0){
                    semiPrimes[i]=1;    
                }
            }
        }
        
        int length = P.length;
        int[] result = new int[length];
        int[] semiprimesAggreation = new int[N+1];
        
        for(int i=1;i<N+1;i++) {
            semiprimesAggreation[i] = semiPrimes[i];
            semiprimesAggreation[i] += semiprimesAggreation[i-1];    
        }
        
        for(int i=0;i<length;i++) { 
            result[i] = semiprimesAggreation[Q[i]] - semiprimesAggreation[P[i]] + semiPrimes[P[i]];
        }
        return result;
    }

    public static int chocolates(int N, int M){
        int gcd = gcd(N, M);
        return N/gcd;
    }

    
    public static int commonPrimeDivs(int[] A, int[] B){
        int cdv = 0;

        for(int i=0; i<A.length; i++){
                if(hasSamePrimeDivisors(A[i], B[i])){
                cdv++;
            }
        }
        return cdv;
    }

    public static boolean hasSamePrimeDivisors(int a, int b) {
	    int gcdValue = gcd(a,b);
        int gcdA;
        int gcdB;
        while(a!=1) {
            gcdA = gcd(a,gcdValue);
            if(gcdA==1)
                break;
            a = a/gcdA;
        }
        if(a!=1)  {
            return false;
        }
        while(b!=1) {
            gcdB = gcd(b,gcdValue);
            if(gcdB==1)
                break;
            b = b/gcdB;
        }
        return b==1;        
	}

    public static int lcm(int a, int b){
        return (a*b)/(gcd(a, b));
    }

    public static int gcd(int a, int b){
        if(a%b==0) return b;
        return gcd(b, a%b);
    }

    public static int fibFrog(int[] A){

        // note: cannot use "List" (both java.util.* and java.awt.* have "List")
        ArrayList<Integer> fibonacci = new ArrayList<>();
        fibonacci.add(0); // note: f(0) = 0 (as in the quesion)
        fibonacci.add(1);
        // note: using "while" is better than "for" (avoid errors)
        while(true){
            int temp1 = fibonacci.get( fibonacci.size()-1 );
            int temp2 = fibonacci.get( fibonacci.size()-2 );
            fibonacci.add( temp1 + temp2 ); 
            
            // if already bigger than length, then break;
            if(temp1 + temp2 > A.length){ 
                break;
            }
        }
        
        // reverse "List": from big to small
        Collections.reverse(fibonacci);
        
        // use "queue" with "point"
        // point(x,y) = point("position", "number of steps")
        ArrayList<Point> queue = new ArrayList<>(); 
        queue.add( new Point(-1, 0) ); // position:-1, steps:0
        
        // index: the current index for queue element 
        int index=0; 
        while(true){
            // cannot take element from queue anymore
            if(index == queue.size() ){ 
                return -1;
            }
            
            // take element from queue
            Point current = queue.get(index); 
            
            // from big to small 
            for(Integer n: fibonacci){
                int nextPosition = current.p + n;
                
                // case 1: "reach the other side" 
                if(nextPosition == A.length){ 
                    // return the number of steps
                    return current.c + 1; 
                }
                
                // case 2: "cannot jump"
                // note: nextPosition < 0 (overflow, be careful)
                else if( (nextPosition > A.length) || (nextPosition < 0)|| (A[nextPosition]==0) ){
                    // note: do nothing 
                }
                
                // case 3: "can jump" (othe cases)
                else{
                    // jump to next position, and step+1
                    Point temp = new Point(nextPosition, current.c + 1); 
                    // add to queue
                    queue.add(temp);  
                    
                    A[nextPosition] = 0; // key point: for high performance~!!
                } 
            }
            
            index++; // take "next element" from queue
        }
    }

    public static int fibonaci(int n){
        if (n<=1){
            return n;
        }
        return fibonaci(n-1) + fibonaci(n - 2);
    }

    public static void main(String argsr[]) {
        int[] a = {15, 10, 3};
        int[] b = {75, 30, 5};

        int[] peeks = {3, 4, 3, 4, 1, 2};

        int[] divs = {3, 1, 2, 3, 6};
        
        String S = "{[()()]}";

        int[] r = {0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0};

        int[] r1 = {0, 1, 0};

        int[] r2 = {1, 1, 1, 0, 1, 0, 0, 0};

        System.out.println(fibFrog(r2));
    }   

}

class Point{
    int p;
    int c;
    public Point(int p, int c){
        this.p=p;
        this.c=c;
    }
}