Pairwise
========


The application programming interface for tcrdist2 involves a set of step-wise
commands centered around the :py:class:`tcrdist.repertoire.TCRrep` class, which
stores input data, methods, user-parameters, and results. It calls methods for
parallelized alignment and distance computation between amino acid strings in the
:py:mod:`tcrdist.pairwise` module.


Although the user does not call any of these functions directly,
this page documents the key methods in the :py:mod:`tcrdist.pairwise` module, upon
which the :py:class:`tcrdist.repertoire.TCRrep` class depends.




Source code: (`pairwise.py <https://github.com/kmayerb/tcrdist2/blob/API2/tcrdist/pairwise.py>`_)

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: tcrdist.pairwise

Internals
---------

.. autofunction:: apply_pw_distance_metric_w_multiprocessing

.. autofunction:: _f_pwdist_parallel_using_distance_wrapper

.. autofunction:: nw_metric

.. autofunction:: hm_metric

.. autofunction:: function_factory

.. autofunction:: _partition

.. autofunction:: get_pwdist_indices

.. autofunction:: get_chunked_pwdist_indices

.. autofunction:: _pack_matrix

.. autofunction:: flatten


Inside Baseball: How is parallel processing used in tcrdist2?
-------------------------------------------------------------

The main function that implements pairwise comparison is
:py:func:`apply_pw_distance_metric_w_multiprocessing`

It take as input a list sequences and a number of control parameters that can be
reviewed in the docstring.


Regardless of the control parameters specified, there are some
common operations as :py:func:`apply_pw_distance_metric_w_multiprocessing`
turns the cranks.

Basically the function does a lot to prepare for one big call:

.. code-block:: python

  multiprocessed_result  = p.map(f, indices)

- **f** is a function that compares two strings from a list of strings.
  by default it will be :py:meth:`_f_pwdist_parallel_using_distance_wrapper`

  .. code-block:: python

    output_tuples = []
      for i,j in indices:
          d = distance_wrapper(unique_seqs[i], unique_seqs[j])
          output_tuples.append((i,j,d, unique_seqs[i], unique_seqs[j]))
      return(output_tuples)


- NOTICE what the **f** is expecting (a list called unique_seqs and a function
  called distance_wrapper in its namespace)

- AND **indicies** a list of lists, each containing of tuples that specify
  which two strings to compare.

- Despite the naming confusiog, indices are mapped to the function **f** across many processes.
  Each **f** will only get one these lists.

But before the fireworks can happen, the function :py:func:`apply_pw_distance_metric_w_multiprocessing`
set's the stage in two ways. It calls a function that creates equal size chunks of
indices and it prepares the namespace of the parallel processes to run :py:meth:`_f_pwdist_parallel_using_distance_wrapper`

This happens in four steps:

#. First, the list of sequences are passed to :py:func:`get_chunked_pwdist_indices`
    - the job of :py:func:`get_chunked_pwdist_indices` is to produce chunks or
      tuples [(i,j),(i,j)] that break up all required comparisons into parts.

#. Next, a distance_wrapper function is created by :py:func:`function_factory` unless
   the heroic user supplies a custom :py:attr:`user_function`. The distance_wrapper is function
   that returns a distance between two strings.  For instance distance_wrapper will most likely
   wrap :py:func:`nw_metric` or :py:func:`hm_metric`.

#. Next, The multiprocessing.Pool() is initialized. Let's look at what goes into this:

      .. code-block:: python

        p = multiprocessing.Pool(processes = processes,
                                 initializer = set_global,
                                 initargs=(sequences,distance_wrapper_new,))



    - set_global() takes initargs

      .. code-block:: python

        def set_global(x,dw):
              global unique_seqs
              unique_seqs = x

              global distance_wrapper
              distance_wrapper = dw

    -  We can see that set_global throws list of sequences into the global namespace
       so each parallel process will have the full list of sequences named as 'unique_seqs',
       which is what :py:func:`_f_pwdist_parallel_using_distance_wrapper`
       is expecting

    -  set_global() also places into the global namespace of the parallel
       a function called 'distance_wrapper' that was
       defined to the user's liking above in the early part of
       :py:func:`apply_pw_distance_metric_w_multiprocessing`


#. All that work leads to the big call:

      .. code-block:: python

        multiprocessed_result  = p.map(f, indices)



Finally because the multiprocessed_result is a list of lists we can return it as a single list
with

  .. code-block:: python

    return(flatten(multiprocessed_result))


At this point :py:func:`apply_pw_distance_metric_w_multiprocessing` is finished,
and the result is sent back to: :py:func:`tcrdist.repertoire._compute_pairwise()`,
which was invoked by :py:class:`tcrdist.repertoire.TCRrep` using the method
:py:meth:`TCRrep.compute_pairwise_all`
