# 15618-Final-Project

website & description: https://tianqi-wu.github.io/course-projects/index.html
Dataset: https://archive.ics.uci.edu/ml/datasets/abalone
This is data of abalone from UCI.

</div>


How to run it: type the following in your console/terminal:

```
make
./knn
```

## description:

<div class="container">
    <h1> Machine Learning Algorithms with Parallelism and Tradeoff: KNN and K-means with OpenMP </h1>
    <h2> Project by Tianqi (Andy) Wu and Alex Zhang</h2>
    <p>URL: tianqi-wu.github.io/course-projects/index.html</p>
    <p>Summary: Our team hopes to explore on how parallelization can improve ML algorithms by building efficient parallelized KNN and K-means clustering with OpenMP, and analyzing tradeoffs based on them.</p>


    <h3>BACKGROUND (from 10-601 lecture notes):</h3>

    <p>Machine learning has been a hot topic in various science and engineering fields. While researchers and engineers admire their accuracy, 
      one can't deny the fact that some algorithms can be really slow to implement when running serially, given the size of data.
      However, such algorithms, when implementing with parallelism, can possibly alleviate such efficiency issues. 
      The authors would like to dive deeper into both supervised and unsupervised algorithms, and implementing one specific algorithm for each one with OpenMP.
    </p>

    <p>For supervised learning, the authors decided to explore K-nearest neighbors algorithm. 
      Such algorithm will label a candidate data point by referring to the K nearest "neighboring" data points. </p>

    <p>For unsupervised learning, the authors decided to explore K-means clustering algorithm. 
      The goal of clustering algorithms in machine learning is to partition unlabeled data into groups with similar features. 
    It can be very useful when analyzing data and classify them. The goal of K-means clustering, in this scenario, is to find 
    an assignment of points to each cluster. </p>

    <p>We would like to explore how parallelism would scale these two algorithms. To be specific, after parallelization is implemented, 
      how the processors may scale these algorithms, etc.</p>
      
    <h3>THE CHALLENGE</h3>
    
    <p>1. The curse of dimensionality. Given the mass amount of data and the possible size of the schema (if any), 
    the data can be relatively difficult to handle, and this can lead to even exponentially growing search time.</p>

    <p>2. Initialization issues. K-means clustering, when not being initialized correctly, can lead to problems. For example, an outlier can 
    drive the result unreasonable from the first glance.</p>

    <p>3. Termination. In the realm of K-means, when and how to terminate the task can be an issue.</p>

    <p>4. Parameters. Whether the parameter is set up correctly may also influence the algorithm.</p>

    <h3>Resources:</h3>

    We hope that we can start the project from scratch. However, we hope that we will refer to existing pesudocode of KNN and K-means clustering,
    
    We also hope that we can use assignment 3 as a reference for general syntax / use cases. We would like to refer to materials such as slides from 10-601,
    
    and possibly papers on what and how to parallelize.

    <h3>Goals</h3>

  <p>1. MUST-DO[50%]: Implement a sequential KNN / K-means clustering algorithms. </p>

  <p>2. MUST-DO[75%]: Implement parallel KNN & K-means algorithms with OpenMP to reduce computational/clustering time. We are thinking about parallelization based on 
    each cluster center / each data point. We hope to scale the code sublinearly at best, given the overhead for assigning and dealing with tasks.</p>

  <p>3. MUST-DO[100%]: record performance on both algorithms with various parameters/hyperparameters and discuss decision tradeoffs.</p>

  <p>4. HOPE-TO-DO[110%]: Find a strategy to aovid early initialization errors (if any) / try various initialization strategies.</p>

  <p>5. HOPE-TO-DO[115%]: Present the visualization process using tools so that it can be easily identifiable.</p>

  <p>6. HOPE-TO-DO[125%]: Find a way to preprocess (if any) to significantly reduce find-nearest time.</p>

    <h3>Platform:</h3>

    We hope to use GHC/PSC clusters. Such platforms provide firm resources of CPU with abundant threads, so that we don't have to worry too much about OpenMP applications.
    Moreover, we have been proficient with such platforms in the past months, and this would give us an advantage.

    <h3>Schedule:</h3>

  <p>4/2 - 4/9: implement KNN and K-means algorithms, sequentially</p>

  <p> 4/9 - 4/16: add OpenMP to naive KNN and K-means, naive performance evaluation</p>

  <p>4/16 - 4/23: write milestone report; advanced performance evaulation (e.g. hyperparameter tuning/initialization), visualization for poster session/final project</p>
    
  <p>4/23 - 4/30: Wrap up and constuct report (possibly accompanying with dynamic visualizations)</p>
    
  <p>4/30 - 5/5: Submit our project, Make a poster with our best effort to make it beautiful and informative.</p>

  <p>The updates for each steps will be posted here accordingly. Thanks for reading this!</p>
