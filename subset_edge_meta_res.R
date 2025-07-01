# the output from step 3
nets = consume_workflow_task_output(wf = 'wf-9654ac7614', 
                                    task_id = NULL,
                                    task_name = 'meta-pvalue-post')

filter_task = cytoreason.cc.client::make_dag_function_task(FUN = function(nets) {
  filtered_output <- nets[nets$cor > 0, ]
  return(filtered_output)
}, 
nets = nets,
task_name ='meta-pvalue-post', 
memory_request = '200Gi',
image =  "eu.gcr.io/cytoreason/ci-cytoreason.individual.variation-package:master_latest"
)

wf <- cytoreason.cc.client::post_dag_workflow(list(filter_task), data_access = 'public')
#  wf-d896d61ba5