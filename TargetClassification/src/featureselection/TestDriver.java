/**
 * 
 */
package featureselection;

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
		TreeGrowthControl ctrl = new TreeGrowthControl();
		ctrl.minCriterion = -Double.MAX_VALUE;
		ctrl.isClassificationUsed = true;
		ctrl.isReplacementUsed = false;
		ctrl.numberOfTreesToGrow = 100;
		int gaRepetitions = 10;
		boolean isXValUsed = true;
		new Controller(args, ctrl, gaRepetitions, isXValUsed);
	}

}
