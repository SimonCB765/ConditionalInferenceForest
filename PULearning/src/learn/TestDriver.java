/**
 * 
 */
package learn;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import tree.ProcessDataForGrowing;
import tree.TreeGrowthControl;

/**
 * @author Simon Bull
 *
 */
public class TestDriver
{

	/**
	 * @param args
	 */
	public static void main(String[] args)
	{

		int testsFailed = 0;
		TreeGrowthControl ctrl = new TreeGrowthControl();
		ctrl.minCriterion = -Double.MAX_VALUE;
		ctrl.isClassificationUsed = true;
		ctrl.isProbabilisticPrediction = true;
		ctrl.numberOfTreesToGrow = 100;
		ctrl.save("C:/Users/Simon/Desktop/PROCDATASAVE.txt");
		TreeGrowthControl ctrlLoad = new TreeGrowthControl("C:/Users/Simon/Desktop/PROCDATASAVE.txt");
		for (Field f : TreeGrowthControl.class.getFields())
		{
			try
			{
				if (!f.get(ctrl).equals(f.get(ctrlLoad)))
				{
					System.out.println(f.getName());
					testsFailed += 1;
				}
			}
			catch (IllegalAccessException e)
			{
				System.err.println("Error: " + e.getMessage());
				System.exit(0);
			}
		}
		ProcessDataForGrowing procData = new ProcessDataForGrowing(args[0], ctrl);
		procData.save("C:/Users/Simon/Desktop/PROCDATASAVE.txt");
		ProcessDataForGrowing procDataLoad = new ProcessDataForGrowing("C:/Users/Simon/Desktop/PROCDATASAVE.txt");
		for (Field f : ProcessDataForGrowing.class.getFields())
		{
			try
			{
				if (!f.get(procData).equals(f.get(procDataLoad)))
				{
					System.out.println(f.getName());
					testsFailed += 1;
				}
			}
			catch (IllegalAccessException e)
			{
				System.err.println("Error: " + e.getMessage());
				System.exit(0);
			}
		}
		System.out.format("Number of tests failed: %d.\n", testsFailed);
		System.exit(0);


		List<List<Integer>> negativesFound = new ArrayList<List<Integer>>();
		for (int i = 0; i < 1; i++)
		{
			System.out.format("This is round %d.\n", i);
			Learner puLearner = new Learner(args[0], ctrl);
			negativesFound.add(puLearner.learnNegativeSet(5, 40));
		}

		Map<Integer, Integer> numNegativeOccurences = new HashMap<Integer, Integer>();
		for (List<Integer> l : negativesFound)
		{
			for (Integer i : l)
			{
				if (numNegativeOccurences.containsKey(i))
				{
					numNegativeOccurences.put(i, numNegativeOccurences.get(i) + 1);
				}
				else
				{
					numNegativeOccurences.put(i, 1);
				}
			}
		}
		System.out.println(numNegativeOccurences.entrySet());

	}

}
