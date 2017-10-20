package viterbi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
public class Viterbi {
	// creates a new Viterbi-type gsm  
	double[] pi;
	double[][] A, B, score;
	int[][] pred;
	int N;
	
	public Viterbi (double[] pi, double[][] A, double[][] B) {
		this.pi = pi;
		this.A = A;
		this.B = B;
		this.N = A.length;
		
		if (errorCheck(pi, A, B)) {
			throw new IllegalArgumentException();
		}
		
		
	} 
	
	boolean errorCheck(double[] pi, double[][] A, double[][] B) {
		if (A[0].length != B.length || A.length < 1 || B[0].length < 1) {
			return true;
		}
		
		if (pi.length != A.length) {
			return true;
		}
		
		double sum = 0;
		for (double d: pi) {
			if (d < 0) {
				return true;
			}
			sum += d;
		}
		if (sum != 1) {
			return true;
		}
		
		//check A, B validity
		if (invalidMatrix(A) || invalidMatrix(B)) {
			return true;
		}
		
		return false;
	}
	
	boolean invalidMatrix(double[][] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			double sum = 0;
			for (int j = 0; j < matrix[0].length; j++) {
				if (matrix[i][j] < 0) {
					return true;
				}
				sum += matrix[i][j];
			}
			if (sum != 1.0) {
				return true;
			}
		}
		return false;
	}
	
	public double probOfSeq(int[] omega, List<Integer> seq, boolean underflow) {
		if (seq == null || seq.size() == 0) {
			throw new IllegalArgumentException();
		}
		
		if (omega == null || omega.length == 0 || omega.length != seq.size()) {
			throw new IllegalArgumentException();
		}
		
		double prob;
		prob = underflow ? Math.log(pi[seq.get(0)]) + 
								Math.log(B[seq.get(0)][omega[0]]) :
						   pi[seq.get(0)] * B[seq.get(0)][omega[0]];
		if (underflow) {
			for (int i = 1; i < seq.size(); i++) {
				prob += Math.log(A[seq.get(i-1)][seq.get(i-1)]) + 
						Math.log(B[seq.get(i)][omega[i]]);
			}
		} else {
			for (int i = 1; i < seq.size(); i++) {
				prob *= A[seq.get(i-1)][seq.get(i-1)] * B[seq.get(i)][omega[i]];
			}	
		}

		return prob;
	}
	
	// computes and prints the sequence with the highest probability of occurring
	public List<Integer> bestSeq(int[] omega, boolean underflow) {
		//forward step of Viterbi algorithm with bestSeq = true
		forwardStep(omega, true, underflow);
		
		//backtrack with bestSeq = true
		LinkedList<Integer> best = backtrack(score, pred, true);
		
		//print results
		printArray(best);
		
		return best;
	}
	
	// computes and prints the sequence with the lowest probability of occurring
	public List<Integer> worstSeq(int[] omega, boolean underflow) {
		// forward step of viterbi algorithm with bestSeq = false
		forwardStep(omega, false, underflow);
		
		//backtrack with bestSeq = false
		LinkedList<Integer> worst = backtrack(score, pred, false);
		
		printArray(worst);
		
		return worst;
	}
	
	void forwardStep(int[] omega, boolean best, boolean underflow) {
		if (omega == null || omega.length == 0) {
			throw new IllegalArgumentException();
		}
		int T = omega.length;
		
		this.score = new double[N][T];
		this.pred = new int[N][T];
		
		for (int i = 0; i < N; i++) {
			score[i][0] = underflow ? Math.log(pi[i]) + Math.log(B[i][omega[0]])
									 : pi[i] * B[i][omega[0]];
		}
		// forward phase
		for (int t = 1; t < T; t++) {
			for (int j = 0; j < N; j++) {
				double tScore = best ? -1 : Double.MAX_VALUE;
				int arg = 0;
				for (int k = 0; k < N; k++) {
					double temp;
					if (underflow) {
						temp = score[k][t-1] + Math.log(A[k][j]) 
											+ Math.log(B[j][omega[t]]);
					} else {
						temp = score[k][t-1] * A[k][j] * B[j][omega[t]];
					}
					if (best) {
						if (temp > tScore) {
							tScore = temp;
							arg = k;
						}
					} else {
						if (temp < tScore) {
							tScore = temp;
							arg = k;
						}
					}
				}
				score[j][t] = tScore;
				pred[j][t] = arg;
			}
		}
	}
	
	LinkedList<Integer> backtrack(double[][] score, int[][] pred, boolean best) {
		LinkedList<Integer> result = new LinkedList<>();
		int T = score[0].length;
		
		int index = 0;
		double firstScore = score[0][T-1];
		
		for (int i = 1; i < N; i++) {
			if (best) {
				index = score[i][T-1] > firstScore ? i : index;
			} else {
				index = score[i][T-1] < firstScore ? i : index;
			}
		}
		
		for (int t = T-1; t >= 1; t--) {
			result.addFirst(index);
			index = pred[index][t];
		}
		result.addFirst(index);
		
		return result;
	}
	
	void printArray(List<Integer> list) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < list.size(); i++) {
			sb.append(list.get(i).toString() + ", ");
		}
		String result = sb.toString();
		System.out.println(result.substring(0, result.lastIndexOf(",")));
	}
	
	

	public List<List<Integer>> bestSeqSet(int[] omega, boolean underflow) {
		int N = this.A.length;
		int T = omega.length;
		
		this.score = new double[N][T];
		int[][] predl = new int[N][T];
		int[][] predr = new int[N][T];
		
		for (int i = 0; i < N; i++) {
			score[i][0] = underflow ? Math.log(pi[i]) + Math.log(B[i][omega[0]])
									 : pi[i] * B[i][omega[0]];
		}
		// forward phase
		for (int t = 1; t < T; t++) {
			for (int j = 0; j < N; j++) {
				double tScore = -1;
				int argMaxL = 0;
				int argMaxR = 0;
				for (int k = 0; k < N; k++) {
					double temp;
					if (underflow) {
						temp = score[k][t-1] + Math.log(A[k][j]) 
											+ Math.log(B[j][omega[t]]);
					} else {
						temp = score[k][t-1] * A[k][j] * B[j][omega[t]];
					}
					if (temp > tScore) {
						tScore = temp;
						argMaxL = k;
						argMaxR = argMaxL;
					} else if (temp == tScore) {
						argMaxR = k;
					}
				}
				score[j][t] = tScore;
				predl[j][t] = argMaxL;
				predr[j][t] = argMaxR;
			}
		}
		List<List<Integer>> result = new ArrayList<>();
		LinkedList<Integer> left = new LinkedList<>();
		LinkedList<Integer> right = new LinkedList<>();
		
		int lIndex = 0;
		int rIndex = 0;
		boolean diff = false;
		double firstScore = score[0][T-1];
		
		for (int i = 1; i < N; i++) {
			if (score[i][T-1] >= firstScore) {
				lIndex = i;
				if (score[i][T-1] == firstScore) {
					rIndex = i;
					diff = true;
				}
			}
		}
		
		if (diff) {
			for (int t = T-1; t >= 1; t--) {
				left.addFirst(lIndex);
				right.addFirst(rIndex);
				lIndex= predl[lIndex][t];
				rIndex = predr[rIndex][t];
			}
			left.addFirst(lIndex);
			right.addFirst(rIndex);
			
			result.add(left);
			result.add(right);
		} else {
			for (int t = T-1; t >= 1; t--) {
				left.addFirst(lIndex);
				lIndex= predl[lIndex][t];
			}
			left.addFirst(lIndex);
			result.add(left);
		}
		
		int cnt = 1;
		for (List<Integer> l: result) {
			if (l != null) {
				System.out.println("List " + cnt);
				printArray(l);
				cnt++;
			}
		}
		
		return result;
	}
	// computes and prints the highest probability found at time T.
	
	public double maxScore(int[] omega, boolean underflow) {
		int T = omega.length;
		forwardStep(omega, true, underflow);
		
		double maxScore = score[0][T-1];
		for (int i = 1; i < N; i++) {
			maxScore = maxScore < score[i][T-1] ? score[i][T-1] : maxScore;
		}
		
		return maxScore;
		
	}
	// computes and prints the k distinct highest probabilities found at time T
	
	public double[] maxScoreSet(int[] omega, int k, boolean underflow) {
		if (omega == null || omega.length == 0 || k < 0) {
			throw new IllegalArgumentException();
		}
		int T = omega.length;
		PriorityQueue<Double> heap = 
				new PriorityQueue<>(k, Collections.reverseOrder());
		
		forwardStep(omega, true, underflow);
		
		for (int i = 0; i < N; i++) {
			heap.add(score[i][T-1]);
		}
		if (k > heap.size()) {
			throw new IllegalArgumentException();
		}
		double[] maxes = new double[k];
		for (int i = 0; i < k; i++) {
			maxes[i] = heap.poll();
		}
		
		System.out.println("Max Score Set of size " + k + ":");
		for (double d: maxes) {
			System.out.println(d);
		}
		System.out.println();
		
		return maxes;
	}
	
	@SuppressWarnings("unchecked")
	public double[] maxScoreSetEC(int[] omega, int k, boolean underflow) {
		if (omega == null || omega.length == 0 || k < 0) {
			throw new IllegalArgumentException();
		}
		int N = this.A.length;
		int T = omega.length;
		
		this.score = new double[N][T];
		
		PriorityQueue<Double>[] pqList = new PriorityQueue[N];
		for (int i = 0; i < N; i++) {
			pqList[i] = new PriorityQueue<>(Collections.reverseOrder());
			pqList[i].add(pi[i] * B[i][omega[0]]);
		}
		
		// forward phase
		for (int t = 1; t < T; t++) {
			PriorityQueue<Double>[] copy = new PriorityQueue[N];
			for (int j = 0; j < N; j++) {
				PriorityQueue<Double> tempPQ = 
						new PriorityQueue<>(Collections.reverseOrder());
				for (int i = 0; i < N; i++) {
					ArrayList<Double> temp = new ArrayList<>();
					while (!pqList[j].isEmpty()) {
						double score = pqList[j].poll();
						double tScore = score * A[i][j] * B[j][omega[t]];
						if (tempPQ.contains(tScore)) {
							tempPQ.add(tScore);
						}
						temp.add(score);
					}
					pqList[i].addAll(temp);
				}
				copy[j] = tempPQ;
			}
			pqList = copy;
		}
		
		PriorityQueue<Double> total = new PriorityQueue<>();
		for (PriorityQueue<Double> pq: pqList) {
			total.addAll(pq);
		}
		
		if (k > total.size()) {
			throw new IllegalArgumentException();
		}
		
		double[] result = new double[k];
		for (int i = 0; i < k; i++) {
			result[i] = total.poll();
		}
		
		for (double d: result) {
			System.out.println(d);
		}
		
		return result;
	}
}