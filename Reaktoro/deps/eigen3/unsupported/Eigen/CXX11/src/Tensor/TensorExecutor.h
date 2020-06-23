// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2014 Benoit Steiner <benoit.steiner.goog@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_CXX11_TENSOR_TENSOR_EXECUTOR_H
#define EIGEN_CXX11_TENSOR_TENSOR_EXECUTOR_H

namespace Eigen {

/**
 * \class TensorExecutor
 * \ingroup CXX11_Tensor_Module
 *
 * \brief The tensor executor class.
 *
 * This class is responsible for launch the evaluation of the expression on
 * the specified computing device.
 *
 * @tparam Vectorizable can use packet math (SSE/AVX/etc... registers and
 *                      instructions)
 * @tparam Tiling       can use block based tensor evaluation
 *                      (see TensorBlock.h)
 */
namespace internal {

/**
 * Evaluating TensorBroadcastingOp via coefficient of packet path is extremely
 * expensive. If expression has at least one broadcast op in it, and it supports
 * block based evaluation, we always prefer it, even for the small tensors. For
 * all other tileable ops, block evaluation overhead for small tensors (fits
 * into L1) is too large, and we fallback on vectorized evaluation.
 */

// TODO(ezhulenev): Add specializations for all other types of Tensor ops.

template<typename Expression>
struct ExpressionHasTensorBroadcastingOp {
  enum { value = false };
};

template<typename LhsXprType, typename RhsXprType>
struct ExpressionHasTensorBroadcastingOp<
    const TensorAssignOp<LhsXprType, RhsXprType> > {
  enum { value = ExpressionHasTensorBroadcastingOp<RhsXprType>::value };
};

template<typename UnaryOp, typename XprType>
struct ExpressionHasTensorBroadcastingOp<
    const TensorCwiseUnaryOp<UnaryOp, XprType> > {
  enum { value = ExpressionHasTensorBroadcastingOp<XprType>::value };
};

template<typename BinaryOp, typename LhsXprType, typename RhsXprType>
struct ExpressionHasTensorBroadcastingOp<
    const TensorCwiseBinaryOp<BinaryOp, LhsXprType, RhsXprType> > {
  enum {
    value = ExpressionHasTensorBroadcastingOp<LhsXprType>::value ||
        ExpressionHasTensorBroadcastingOp<RhsXprType>::value
  };
};

template<typename Broadcast, typename XprType>
struct ExpressionHasTensorBroadcastingOp<
    const TensorBroadcastingOp<Broadcast, XprType> > {
  enum { value = true };
};

// -------------------------------------------------------------------------- //

/**
 * Default strategy: the expression is evaluated sequentially with a single cpu
 * thread, without vectorization and block evaluation.
 */
<<<<<<< HEAD
template <typename Expression, typename Device, bool Vectorizable,
          TiledEvaluation Tiling>
=======
#if EIGEN_HAS_CXX11
template <typename Expression, typename Device, bool Vectorizable,
          TiledEvaluation Tiling>
#else
 template <typename Expression, typename Device, bool Vectorizable,
          TiledEvaluation::TiledEvaluation Tiling>
#endif
>>>>>>> master
class TensorExecutor {
 public:
  typedef typename Expression::Index StorageIndex;

<<<<<<< HEAD
  // Including `unsupported/Eigen/CXX11/Tensor` in different translation units
  // with/without `EIGEN_USE_THREADS` or `EIGEN_USE_GPU` is a potential ODR
  // violation. If this template is instantiated with a non-default device, it
  // means that this header file was included without defining
  // `EIGEN_USE_THREADS`, `EIGEN_USE_GPU` or `EIGEN_USE_SYCL`.
  static_assert(std::is_same<Device, DefaultDevice>::value,
                "Default executor instantiated with non-default device. "
                "You must #define EIGEN_USE_THREADS, EIGEN_USE_GPU or "
                "EIGEN_USE_SYCL before including Eigen headers.");

=======
>>>>>>> master
  EIGEN_DEVICE_FUNC
  static EIGEN_STRONG_INLINE void run(const Expression& expr,
                                      const Device& device = Device()) {
    TensorEvaluator<Expression, Device> evaluator(expr, device);
    const bool needs_assign = evaluator.evalSubExprsIfNeeded(NULL);
    if (needs_assign) {
      const StorageIndex size = array_prod(evaluator.dimensions());
      for (StorageIndex i = 0; i < size; ++i) {
        evaluator.evalScalar(i);
      }
    }
    evaluator.cleanup();
  }
};

/**
 * Default async execution strategy is not implemented. Currently it's only
 * available for ThreadPoolDevice (see definition below).
 */
template <typename Expression, typename Device, typename DoneCallback,
<<<<<<< HEAD
          bool Vectorizable, TiledEvaluation Tiling>
=======
          bool Vectorizable, bool Tileable>
>>>>>>> master
class TensorAsyncExecutor {};

/**
 * Process all the data with a single cpu thread, using vectorized instructions.
 */
template <typename Expression>
class TensorExecutor<Expression, DefaultDevice, /*Vectorizable=*/true,
                     /*Tiling=*/TiledEvaluation::Off> {
 public:
  typedef typename Expression::Index StorageIndex;

  EIGEN_DEVICE_FUNC
  static EIGEN_STRONG_INLINE void run(
      const Expression& expr, const DefaultDevice& device = DefaultDevice()) {
    TensorEvaluator<Expression, DefaultDevice> evaluator(expr, device);
    const bool needs_assign = evaluator.evalSubExprsIfNeeded(NULL);
    if (needs_assign) {
      const StorageIndex size = array_prod(evaluator.dimensions());
      const int PacketSize = unpacket_traits<typename TensorEvaluator<
          Expression, DefaultDevice>::PacketReturnType>::size;

      // Give compiler a strong possibility to unroll the loop. But don't insist
      // on unrolling, because if the function is expensive compiler should not
      // unroll the loop at the expense of inlining.
      const StorageIndex UnrolledSize =
          (size / (4 * PacketSize)) * 4 * PacketSize;
      for (StorageIndex i = 0; i < UnrolledSize; i += 4 * PacketSize) {
        for (StorageIndex j = 0; j < 4; j++) {
          evaluator.evalPacket(i + j * PacketSize);
        }
      }
      const StorageIndex VectorizedSize = (size / PacketSize) * PacketSize;
      for (StorageIndex i = UnrolledSize; i < VectorizedSize; i += PacketSize) {
        evaluator.evalPacket(i);
      }
      for (StorageIndex i = VectorizedSize; i < size; ++i) {
        evaluator.evalScalar(i);
      }
    }
    evaluator.cleanup();
  }
};

/**
 * Process all the data with a single cpu thread, using blocks of data. By
 * sizing a block to fit L1 cache we get better cache performance.
 */
template <typename Expression, bool Vectorizable>
class TensorExecutor<Expression, DefaultDevice, Vectorizable,
<<<<<<< HEAD
=======
                     /*Tiling=*/TiledEvaluation::Legacy> {
 public:
  typedef typename traits<Expression>::Scalar Scalar;
  typedef typename remove_const<Scalar>::type ScalarNoConst;

  typedef TensorEvaluator<Expression, DefaultDevice> Evaluator;
  typedef typename traits<Expression>::Index StorageIndex;

  static const int NumDims = traits<Expression>::NumDimensions;

  EIGEN_DEVICE_FUNC
  static EIGEN_STRONG_INLINE void run(const Expression& expr,
                         const DefaultDevice& device = DefaultDevice()) {
    typedef TensorBlock<ScalarNoConst, StorageIndex, NumDims, Evaluator::Layout> TensorBlock;
    typedef TensorBlockMapper<ScalarNoConst, StorageIndex, NumDims, Evaluator::Layout> TensorBlockMapper;
    typedef typename TensorBlock::Dimensions TensorBlockDimensions;

    Evaluator evaluator(expr, device);
    Index total_size = array_prod(evaluator.dimensions());
    Index cache_size = device.firstLevelCacheSize() / sizeof(Scalar);

    if (total_size < cache_size
        && !ExpressionHasTensorBroadcastingOp<Expression>::value) {
      // TODO(andydavis) Reduce block management overhead for small tensors.
      internal::TensorExecutor<Expression, DefaultDevice, Vectorizable, /*Tiling=*/TiledEvaluation::Off>::run(expr,device);
      evaluator.cleanup();
      return;
    }

    const bool needs_assign = evaluator.evalSubExprsIfNeeded(NULL);
    if (needs_assign) {
      // Size tensor blocks to fit in cache (or requested target block size).
      Index block_total_size = numext::mini(cache_size, total_size);
      TensorBlockShapeType block_shape = kSkewedInnerDims;
      // Query expression tree for desired block size/shape.
      std::vector<TensorOpResourceRequirements> resources;
      evaluator.getResourceRequirements(&resources);
      MergeResourceRequirements(resources, &block_shape, &block_total_size);

      TensorBlockMapper block_mapper(
          TensorBlockDimensions(evaluator.dimensions()), block_shape,
          block_total_size);
      block_total_size = block_mapper.block_dims_total_size();

      ScalarNoConst* data = static_cast<ScalarNoConst*>(
          device.allocate(block_total_size * sizeof(Scalar)));

      const StorageIndex total_block_count = block_mapper.total_block_count();
      for (StorageIndex i = 0; i < total_block_count; ++i) {
        TensorBlock block = block_mapper.GetBlockForIndex(i, data);
        evaluator.evalBlock(&block);
      }
      device.deallocate(data);
    }
    evaluator.cleanup();
  }
};

/**
 * Process all the data with a single cpu thread, using blocks of data. By
 * sizing a block to fit L1 cache we get better cache performance.
 */
template <typename Expression, bool Vectorizable>
class TensorExecutor<Expression, DefaultDevice, Vectorizable,
>>>>>>> master
                     /*Tiling=*/TiledEvaluation::On> {
 public:
  typedef typename traits<Expression>::Scalar Scalar;
  typedef typename remove_const<Scalar>::type ScalarNoConst;

  typedef TensorEvaluator<Expression, DefaultDevice> Evaluator;
  typedef typename traits<Expression>::Index StorageIndex;

  static const int NumDims = traits<Expression>::NumDimensions;

  EIGEN_DEVICE_FUNC
  static EIGEN_STRONG_INLINE void run(const Expression& expr,
                         const DefaultDevice& device = DefaultDevice()) {
<<<<<<< HEAD
    typedef TensorBlockMapper<NumDims, Evaluator::Layout, StorageIndex>
        TensorBlockMapper;
=======
    typedef TensorBlock<ScalarNoConst, StorageIndex, NumDims, Evaluator::Layout> TensorBlock;
    typedef TensorBlockMapper<ScalarNoConst, StorageIndex, NumDims, Evaluator::Layout> TensorBlockMapper;
    typedef typename TensorBlock::Dimensions TensorBlockDimensions;
>>>>>>> master

    typedef internal::TensorBlockDescriptor<NumDims, StorageIndex>
        TensorBlockDesc;
    typedef internal::TensorBlockScratchAllocator<DefaultDevice>
        TensorBlockScratch;

    Evaluator evaluator(expr, device);
<<<<<<< HEAD
=======
    Index total_size = array_prod(evaluator.dimensions());
    Index cache_size = device.firstLevelCacheSize() / sizeof(Scalar);
>>>>>>> master

    // TODO(ezhulenev): Do not use tiling for small tensors?
    const bool needs_assign = evaluator.evalSubExprsIfNeeded(NULL);

    if (needs_assign) {
<<<<<<< HEAD
      // Query expression tree for desired block size/shape.
      const TensorBlockResourceRequirements requirements =
          evaluator.getResourceRequirements();

      const TensorBlockMapper block_mapper(
          typename TensorBlockDesc::Dimensions(evaluator.dimensions()),
          requirements);
=======
      // Size tensor blocks to fit in cache (or requested target block size).
      Index block_total_size = numext::mini(cache_size, total_size);
      TensorBlockShapeType block_shape = kSkewedInnerDims;
      // Query expression tree for desired block size/shape.
      std::vector<TensorOpResourceRequirements> resources;
      evaluator.getResourceRequirements(&resources);
      MergeResourceRequirements(resources, &block_shape, &block_total_size);

      TensorBlockMapper block_mapper(
          TensorBlockDimensions(evaluator.dimensions()), block_shape,
          block_total_size);
      block_total_size = block_mapper.block_dims_total_size();
>>>>>>> master

      // Share scratch memory allocator between all blocks.
      TensorBlockScratch scratch(device);

<<<<<<< HEAD
      const StorageIndex total_block_count = block_mapper.blockCount();
      for (StorageIndex i = 0; i < total_block_count; ++i) {
        TensorBlockDesc desc = block_mapper.blockDescriptor(i);
        evaluator.evalBlock(desc, scratch);
=======
      const StorageIndex total_block_count = block_mapper.total_block_count();
      for (StorageIndex i = 0; i < total_block_count; ++i) {
        TensorBlock block = block_mapper.GetBlockForIndex(i, NULL);

        TensorBlockDesc desc(block.first_coeff_index(), block.block_sizes());
        evaluator.evalBlockV2(desc, scratch);
>>>>>>> master
        scratch.reset();
      }
    }
    evaluator.cleanup();
  }
};

/**
 * Multicore strategy: the index space is partitioned and each partition is
 * executed on a single core.
 *
 * (1) TensorExecutor will submit work to the ThreadPoolDevice managed thread
 *     pool, and will block the caller thread until all tasks are finished.
 *
 * (2) TensorAsyncExecutor is a non-blocking version, that will submit work to
 *     the ThreadPoolDevice managed thread pool, and will return immediately.
 *     It will call 'done' callback after all tasks are finished.
 */
#ifdef EIGEN_USE_THREADS

template <typename TensorBlockMapper>
struct TensorExecutorTilingContext {
<<<<<<< HEAD
  TensorExecutorTilingContext() = default;
  TensorExecutorTilingContext(const TensorBlockMapper& b_mapper,
                              const TensorOpCost& b_cost, size_t b_aligned_size)
      : block_mapper(b_mapper),
        cost(b_cost),
        aligned_blocksize(b_aligned_size) {}

  TensorBlockMapper block_mapper;  // navigate through blocks
  TensorOpCost cost;               // cost of computing a single block
=======
  typedef typename TensorBlockMapper::Block TensorBlock;

  TensorExecutorTilingContext() : buffer(nullptr) {}
  TensorExecutorTilingContext(const TensorBlockMapper& b_mapper,
                              const TensorOpCost& b_cost, void* b_buffer,
                              size_t b_aligned_size)
      : block_mapper(b_mapper),
        cost(b_cost),
        buffer(b_buffer),
        aligned_blocksize(b_aligned_size) {}

  template <typename Scalar>
  Scalar* GetCurrentThreadBuffer(const ThreadPoolDevice& device) const {
    // ThreadPoolDevice::currentThreadId() returns -1 if called from a thread
    // not in the thread pool, such as the main thread dispatching Eigen
    // expressions.
    const int thread_idx = device.currentThreadId();
    eigen_assert(thread_idx >= -1 && thread_idx < device.numThreads());

    const Index offset = aligned_blocksize * (thread_idx + 1);
    return reinterpret_cast<Scalar*>(static_cast<char*>(buffer) + offset);
  }

  TensorBlockMapper block_mapper;  // navigate through blocks
  TensorOpCost cost;               // cost of computing a single block
  void* buffer;                    // temporary buffer for blocks
>>>>>>> master
  size_t aligned_blocksize;        // block size after memory alignment
};

// Computes a block evaluation parameters, and allocates temporary memory buffer
// for blocks. See TensorExecutor/TensorAsyncExecutor (Tiling=On) below.
template <typename Evaluator, typename TensorBlockMapper, bool Vectorizable>
TensorExecutorTilingContext<TensorBlockMapper> GetTensorExecutorTilingContext(
<<<<<<< HEAD
    const Evaluator& evaluator) {
  // Query expression tree for desired block size/shape.
  TensorBlockResourceRequirements requirements =
      evaluator.getResourceRequirements();

  // Update target block size based on cost model.
  double taskSize = TensorCostModel<ThreadPoolDevice>::taskSize(
      1, requirements.cost_per_coeff);
  requirements.size = static_cast<size_t>(1.0 / taskSize);

  TensorBlockMapper block_mapper(
      typename TensorBlockMapper::Dimensions(evaluator.dimensions()),
      requirements);

  size_t block_size = block_mapper.blockTotalSize();
=======
    const ThreadPoolDevice& device, const Evaluator& evaluator,
    bool allocate_buffer = true) {
  // Prefer blocks skewed toward inner dimension.
  TensorBlockShapeType block_shape = kSkewedInnerDims;
  Index block_total_size = 0;

  // Query expression tree for desired block size/shape.
  std::vector<TensorOpResourceRequirements> resources;
  evaluator.getResourceRequirements(&resources);
  MergeResourceRequirements(resources, &block_shape, &block_total_size);
  int num_threads = device.numThreads();

  // Estimate minimum block size based on cost.
  TensorOpCost cost = evaluator.costPerCoeff(Vectorizable);
  double taskSize = TensorCostModel<ThreadPoolDevice>::taskSize(1, cost);
  size_t block_size = static_cast<size_t>(1.0 / taskSize);

  TensorBlockMapper block_mapper(
      typename TensorBlockMapper::Dimensions(evaluator.dimensions()),
      block_shape, block_size);

  block_size = block_mapper.block_dims_total_size();
>>>>>>> master
  const size_t align = numext::maxi(EIGEN_MAX_ALIGN_BYTES, 1);
  const size_t aligned_blocksize =
      align *
      divup<size_t>(block_size * sizeof(typename Evaluator::Scalar), align);

<<<<<<< HEAD
  return {block_mapper, requirements.cost_per_coeff * block_size,
          aligned_blocksize};
=======
  // TODO(ezhulenev): In new block evaluation framework there is no need for
  // allocating temporary buffers, remove this after migration.
  void* buf = NULL;
  if (allocate_buffer) {
    buf = device.allocate((num_threads + 1) * aligned_blocksize);
  }

  return {block_mapper, cost * block_size, buf, aligned_blocksize};
>>>>>>> master
}

template <typename Evaluator, typename StorageIndex, bool Vectorizable>
struct EvalRange {
  static void run(Evaluator* evaluator_in, const StorageIndex firstIdx,
                  const StorageIndex lastIdx) {
    Evaluator evaluator = *evaluator_in;
    eigen_assert(lastIdx >= firstIdx);
    for (StorageIndex i = firstIdx; i < lastIdx; ++i) {
      evaluator.evalScalar(i);
    }
  }

  static StorageIndex alignBlockSize(StorageIndex size) { return size; }
};

template <typename Evaluator, typename StorageIndex>
struct EvalRange<Evaluator, StorageIndex, /*Vectorizable*/ true> {
  static const int PacketSize =
      unpacket_traits<typename Evaluator::PacketReturnType>::size;

  static void run(Evaluator* evaluator_in, const StorageIndex firstIdx,
                  const StorageIndex lastIdx) {
    Evaluator evaluator = *evaluator_in;
    eigen_assert(lastIdx >= firstIdx);
    StorageIndex i = firstIdx;
    if (lastIdx - firstIdx >= PacketSize) {
      eigen_assert(firstIdx % PacketSize == 0);
      StorageIndex last_chunk_offset = lastIdx - 4 * PacketSize;
      // Give compiler a strong possibility to unroll the loop. But don't insist
      // on unrolling, because if the function is expensive compiler should not
      // unroll the loop at the expense of inlining.
      for (; i <= last_chunk_offset; i += 4 * PacketSize) {
        for (StorageIndex j = 0; j < 4; j++) {
          evaluator.evalPacket(i + j * PacketSize);
        }
      }
      last_chunk_offset = lastIdx - PacketSize;
      for (; i <= last_chunk_offset; i += PacketSize) {
        evaluator.evalPacket(i);
      }
    }
    for (; i < lastIdx; ++i) {
      evaluator.evalScalar(i);
    }
  }

  static StorageIndex alignBlockSize(StorageIndex size) {
    // Align block size to packet size and account for unrolling in run above.
    if (size >= 16 * PacketSize) {
      return (size + 4 * PacketSize - 1) & ~(4 * PacketSize - 1);
    }
    // Aligning to 4 * PacketSize would increase block size by more than 25%.
    return (size + PacketSize - 1) & ~(PacketSize - 1);
  }
};

template <typename Expression, bool Vectorizable, TiledEvaluation Tiling>
class TensorExecutor<Expression, ThreadPoolDevice, Vectorizable, Tiling> {
 public:
  typedef typename Expression::Index StorageIndex;

  static EIGEN_STRONG_INLINE void run(const Expression& expr,
                         const ThreadPoolDevice& device) {
    typedef TensorEvaluator<Expression, ThreadPoolDevice> Evaluator;
    typedef EvalRange<Evaluator, StorageIndex, Vectorizable> EvalRange;

    Evaluator evaluator(expr, device);
    const bool needs_assign = evaluator.evalSubExprsIfNeeded(nullptr);
    if (needs_assign) {
      const StorageIndex size = array_prod(evaluator.dimensions());
      device.parallelFor(size, evaluator.costPerCoeff(Vectorizable),
                         EvalRange::alignBlockSize,
                         [&evaluator](StorageIndex firstIdx, StorageIndex lastIdx) {
                           EvalRange::run(&evaluator, firstIdx, lastIdx);
                         });
    }
    evaluator.cleanup();
  }
};

template <typename Expression, bool Vectorizable>
class TensorExecutor<Expression, ThreadPoolDevice, Vectorizable,
<<<<<<< HEAD
=======
                     /*Tiling=*/TiledEvaluation::Legacy> {
 public:
  typedef typename traits<Expression>::Index StorageIndex;
  typedef typename traits<Expression>::Scalar Scalar;
  typedef typename remove_const<Scalar>::type ScalarNoConst;

  static const int NumDims = traits<Expression>::NumDimensions;

  typedef TensorEvaluator<Expression, ThreadPoolDevice> Evaluator;
  typedef TensorBlockMapper<ScalarNoConst, StorageIndex, NumDims, Evaluator::Layout> BlockMapper;
  typedef TensorExecutorTilingContext<BlockMapper> TilingContext;

  static EIGEN_STRONG_INLINE void run(const Expression& expr,
                         const ThreadPoolDevice& device) {
    Evaluator evaluator(expr, device);
    Index total_size = array_prod(evaluator.dimensions());
    Index cache_size = device.firstLevelCacheSize() / sizeof(Scalar);

    if (total_size < cache_size &&
        !ExpressionHasTensorBroadcastingOp<Expression>::value) {
      // TODO(andydavis) Reduce block management overhead for small tensors.
      internal::TensorExecutor<Expression, ThreadPoolDevice, Vectorizable,
                               /*Tiling=*/TiledEvaluation::Off>::run(expr,
                                                                     device);
      evaluator.cleanup();
      return;
    }

    const bool needs_assign = evaluator.evalSubExprsIfNeeded(nullptr);
    if (needs_assign) {
      const TilingContext tiling =
          internal::GetTensorExecutorTilingContext<Evaluator, BlockMapper,
                                                   Vectorizable>(device, evaluator);

      device.parallelFor(
          tiling.block_mapper.total_block_count(), tiling.cost,
          [=, &device, &evaluator, &tiling](StorageIndex firstIdx,
                                            StorageIndex lastIdx) {
            ScalarNoConst* thread_buf =
                tiling.template GetCurrentThreadBuffer<ScalarNoConst>(device);
            for (StorageIndex i = firstIdx; i < lastIdx; ++i) {
              auto block = tiling.block_mapper.GetBlockForIndex(i, thread_buf);
              evaluator.evalBlock(&block);
            }
          });
      device.deallocate(tiling.buffer);
    }
    evaluator.cleanup();
  }
};

template <typename Expression, bool Vectorizable>
class TensorExecutor<Expression, ThreadPoolDevice, Vectorizable,
>>>>>>> master
                     /*Tiling=*/TiledEvaluation::On> {
 public:
  typedef typename traits<Expression>::Index IndexType;
  typedef typename traits<Expression>::Scalar Scalar;
  typedef typename remove_const<Scalar>::type ScalarNoConst;

  static const int NumDims = traits<Expression>::NumDimensions;

  typedef TensorEvaluator<Expression, ThreadPoolDevice> Evaluator;
<<<<<<< HEAD
  typedef TensorBlockMapper<NumDims, Evaluator::Layout, IndexType> BlockMapper;
=======
  typedef TensorBlockMapper<ScalarNoConst, IndexType, NumDims,
                            Evaluator::Layout>
      BlockMapper;
>>>>>>> master
  typedef TensorExecutorTilingContext<BlockMapper> TilingContext;

  typedef internal::TensorBlockDescriptor<NumDims, IndexType>
      TensorBlockDesc;
  typedef internal::TensorBlockScratchAllocator<ThreadPoolDevice>
      TensorBlockScratch;

  static EIGEN_STRONG_INLINE void run(const Expression& expr,
                                      const ThreadPoolDevice& device) {
    Evaluator evaluator(expr, device);

    const bool needs_assign = evaluator.evalSubExprsIfNeeded(nullptr);
    if (needs_assign) {
      const TilingContext tiling =
          internal::GetTensorExecutorTilingContext<Evaluator, BlockMapper,
<<<<<<< HEAD
                                                   Vectorizable>(evaluator);
=======
                                                   Vectorizable>(
              device, evaluator, /*allocate_buffer=*/false);
>>>>>>> master

      auto eval_block = [&device, &evaluator, &tiling](IndexType firstBlockIdx,
                                                       IndexType lastBlockIdx) {
        TensorBlockScratch scratch(device);

<<<<<<< HEAD
        for (IndexType block_idx = firstBlockIdx; block_idx < lastBlockIdx;
             ++block_idx) {
          TensorBlockDesc desc = tiling.block_mapper.blockDescriptor(block_idx);
          evaluator.evalBlock(desc, scratch);
=======
        for (IndexType block_idx = firstBlockIdx; block_idx < lastBlockIdx; ++block_idx) {
          auto block = tiling.block_mapper.GetBlockForIndex(block_idx, nullptr);
          TensorBlockDesc desc(block.first_coeff_index(), block.block_sizes());
          evaluator.evalBlockV2(desc, scratch);
>>>>>>> master
          scratch.reset();
        }
      };

<<<<<<< HEAD
      // Evaluate small expressions directly as a single block.
      if (tiling.block_mapper.blockCount() == 1) {
        TensorBlockScratch scratch(device);
        TensorBlockDesc desc(0, tiling.block_mapper.blockDimensions());
        evaluator.evalBlock(desc, scratch);
      } else {
        device.parallelFor(tiling.block_mapper.blockCount(), tiling.cost,
                           eval_block);
      }
=======
      device.parallelFor(tiling.block_mapper.total_block_count(), tiling.cost,
                         eval_block);
>>>>>>> master
    }
    evaluator.cleanup();
  }
};

template <typename Expression, typename DoneCallback, bool Vectorizable,
<<<<<<< HEAD
          TiledEvaluation Tiling>
class TensorAsyncExecutor<Expression, ThreadPoolDevice, DoneCallback,
                          Vectorizable, Tiling> {
=======
          bool Tileable>
class TensorAsyncExecutor<Expression, ThreadPoolDevice, DoneCallback,
                          Vectorizable, Tileable> {
>>>>>>> master
 public:
  typedef typename Expression::Index StorageIndex;
  typedef TensorEvaluator<Expression, ThreadPoolDevice> Evaluator;

  static EIGEN_STRONG_INLINE void runAsync(const Expression& expr,
                                           const ThreadPoolDevice& device,
                                           DoneCallback done) {
    TensorAsyncExecutorContext* const ctx =
        new TensorAsyncExecutorContext(expr, device, std::move(done));

    const auto on_eval_subexprs = [ctx, &device](bool need_assign) -> void {
      if (!need_assign) {
        delete ctx;
        return;
      }

      typedef EvalRange<Evaluator, StorageIndex, Vectorizable> EvalRange;
      const StorageIndex size = array_prod(ctx->evaluator.dimensions());
      device.parallelForAsync(
          size, ctx->evaluator.costPerCoeff(Vectorizable),
          EvalRange::alignBlockSize,
          [ctx](StorageIndex firstIdx, StorageIndex lastIdx) {
            EvalRange::run(&ctx->evaluator, firstIdx, lastIdx);
          },
          [ctx]() { delete ctx; });
    };

    ctx->evaluator.evalSubExprsIfNeededAsync(nullptr, on_eval_subexprs);
  }

 private:
  struct TensorAsyncExecutorContext {
    TensorAsyncExecutorContext(const Expression& expr,
                               const ThreadPoolDevice& thread_pool,
                               DoneCallback done)
        : evaluator(expr, thread_pool), on_done(std::move(done)) {}

    ~TensorAsyncExecutorContext() {
<<<<<<< HEAD
      evaluator.cleanup();
      on_done();
=======
      on_done();
      evaluator.cleanup();
>>>>>>> master
    }

    Evaluator evaluator;

   private:
    DoneCallback on_done;
  };
};

template <typename Expression, typename DoneCallback, bool Vectorizable>
class TensorAsyncExecutor<Expression, ThreadPoolDevice, DoneCallback,
<<<<<<< HEAD
                          Vectorizable, /*Tileable*/ TiledEvaluation::On> {
 public:
  typedef typename traits<Expression>::Index IndexType;
=======
                          Vectorizable, /*Tileable*/ true> {
 public:
  typedef typename traits<Expression>::Index StorageIndex;
>>>>>>> master
  typedef typename traits<Expression>::Scalar Scalar;
  typedef typename remove_const<Scalar>::type ScalarNoConst;

  static const int NumDims = traits<Expression>::NumDimensions;

  typedef TensorEvaluator<Expression, ThreadPoolDevice> Evaluator;
<<<<<<< HEAD
  typedef TensorBlockMapper<NumDims, Evaluator::Layout, IndexType> BlockMapper;
  typedef TensorExecutorTilingContext<BlockMapper> TilingContext;

  typedef internal::TensorBlockDescriptor<NumDims, IndexType> TensorBlockDesc;
  typedef internal::TensorBlockScratchAllocator<ThreadPoolDevice>
      TensorBlockScratch;

  static EIGEN_STRONG_INLINE void runAsync(const Expression& expr,
                                           const ThreadPoolDevice& device,
                                           DoneCallback done) {

    TensorAsyncExecutorContext* const ctx =
        new TensorAsyncExecutorContext(expr, device, std::move(done));

    const auto on_eval_subexprs = [ctx](bool need_assign) -> void {
=======
  typedef TensorBlockMapper<ScalarNoConst, StorageIndex, NumDims, Evaluator::Layout> BlockMapper;
  typedef TensorExecutorTilingContext<BlockMapper> TilingContext;

  static EIGEN_STRONG_INLINE void runAsync(const Expression& expr,
                                           const ThreadPoolDevice& device,
                                           DoneCallback done) {
    TensorAsyncExecutorContext* const ctx =
        new TensorAsyncExecutorContext(expr, device, std::move(done));

    Index total_size = array_prod(ctx->evaluator.dimensions());
    Index cache_size = device.firstLevelCacheSize() / sizeof(Scalar);

    if (total_size < cache_size &&
        !ExpressionHasTensorBroadcastingOp<Expression>::value) {
      auto delete_ctx = [ctx]() { delete ctx; };
      internal::TensorAsyncExecutor<
          Expression, ThreadPoolDevice, decltype(delete_ctx), Vectorizable,
          /*Tileable*/ false>::runAsync(expr, device, std::move(delete_ctx));
      return;
    }

    const auto on_eval_subexprs = [ctx, &device](bool need_assign) -> void {
>>>>>>> master
      if (!need_assign) {
        delete ctx;
        return;
      }

<<<<<<< HEAD
      ctx->tiling = internal::GetTensorExecutorTilingContext<
          Evaluator, BlockMapper, Vectorizable>(ctx->evaluator);

      auto eval_block = [ctx](IndexType firstBlockIdx, IndexType lastBlockIdx) {
        TensorBlockScratch scratch(ctx->device);

        for (IndexType block_idx = firstBlockIdx; block_idx < lastBlockIdx;
             ++block_idx) {
          TensorBlockDesc desc =
              ctx->tiling.block_mapper.blockDescriptor(block_idx);
          ctx->evaluator.evalBlock(desc, scratch);
          scratch.reset();
        }
      };

      // Evaluate small expressions directly as a single block.
      if (ctx->tiling.block_mapper.blockCount() == 1) {
        TensorBlockScratch scratch(ctx->device);
        TensorBlockDesc desc(0, ctx->tiling.block_mapper.blockDimensions());
        ctx->evaluator.evalBlock(desc, scratch);
        delete ctx;
      } else {
        ctx->device.parallelForAsync(ctx->tiling.block_mapper.blockCount(),
                                     ctx->tiling.cost, eval_block,
                                     [ctx]() { delete ctx; });
      }
=======
      ctx->tiling =
          GetTensorExecutorTilingContext<Evaluator, BlockMapper,
                                         Vectorizable>(device, ctx->evaluator);

      device.parallelForAsync(
          ctx->tiling.block_mapper.total_block_count(), ctx->tiling.cost,
          [ctx](StorageIndex firstIdx, StorageIndex lastIdx) {
            ScalarNoConst* thread_buf =
                ctx->tiling.template GetCurrentThreadBuffer<ScalarNoConst>(
                    ctx->device);
            for (StorageIndex i = firstIdx; i < lastIdx; ++i) {
              auto block =
                  ctx->tiling.block_mapper.GetBlockForIndex(i, thread_buf);
              ctx->evaluator.evalBlock(&block);
            }
          },
          [ctx]() { delete ctx; });
>>>>>>> master
    };

    ctx->evaluator.evalSubExprsIfNeededAsync(nullptr, on_eval_subexprs);
  }

 private:
  struct TensorAsyncExecutorContext {
    TensorAsyncExecutorContext(const Expression& expr,
                               const ThreadPoolDevice& thread_pool,
                               DoneCallback done)
        : device(thread_pool),
          evaluator(expr, thread_pool),
          on_done(std::move(done)) {}

    ~TensorAsyncExecutorContext() {
<<<<<<< HEAD
      evaluator.cleanup();
      on_done();
=======
      on_done();
      device.deallocate(tiling.buffer);
      evaluator.cleanup();
>>>>>>> master
    }

    const ThreadPoolDevice& device;
    Evaluator evaluator;
    TilingContext tiling;

   private:
    DoneCallback on_done;
  };
};

#endif  // EIGEN_USE_THREADS

<<<<<<< HEAD
=======

>>>>>>> master
// GPU: the evaluation of the expression is offloaded to a GPU.
#if defined(EIGEN_USE_GPU)

template <typename Expression, bool Vectorizable, TiledEvaluation Tiling>
class TensorExecutor<Expression, GpuDevice, Vectorizable, Tiling> {
 public:
  typedef typename Expression::Index StorageIndex;
  static void run(const Expression& expr, const GpuDevice& device);
};

#if defined(EIGEN_GPUCC)
template <typename Evaluator, typename StorageIndex, bool Vectorizable>
struct EigenMetaKernelEval {
<<<<<<< HEAD
  static EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE
=======
  static __device__ EIGEN_ALWAYS_INLINE
>>>>>>> master
  void run(Evaluator& eval, StorageIndex firstIdx, StorageIndex lastIdx, StorageIndex step_size) {
    for (StorageIndex i = firstIdx; i < lastIdx; i += step_size) {
      eval.evalScalar(i);
    }
  }
};

template <typename Evaluator, typename StorageIndex>
struct EigenMetaKernelEval<Evaluator, StorageIndex, true> {
<<<<<<< HEAD
  static EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE
=======
  static __device__ EIGEN_ALWAYS_INLINE
>>>>>>> master
  void run(Evaluator& eval, StorageIndex firstIdx, StorageIndex lastIdx, StorageIndex step_size) {
    const StorageIndex PacketSize = unpacket_traits<typename Evaluator::PacketReturnType>::size;
    const StorageIndex vectorized_size = (lastIdx / PacketSize) * PacketSize;
    const StorageIndex vectorized_step_size = step_size * PacketSize;

    // Use the vector path
    for (StorageIndex i = firstIdx * PacketSize; i < vectorized_size;
         i += vectorized_step_size) {
      eval.evalPacket(i);
    }
    for (StorageIndex i = vectorized_size + firstIdx; i < lastIdx; i += step_size) {
      eval.evalScalar(i);
    }
  }
};

template <typename Evaluator, typename StorageIndex>
__global__ void
__launch_bounds__(1024)
EigenMetaKernel(Evaluator eval, StorageIndex size) {

  const StorageIndex first_index = blockIdx.x * blockDim.x + threadIdx.x;
  const StorageIndex step_size = blockDim.x * gridDim.x;

  const bool vectorizable = Evaluator::PacketAccess & Evaluator::IsAligned;
  EigenMetaKernelEval<Evaluator, StorageIndex, vectorizable>::run(eval, first_index, size, step_size);
}

/*static*/
template <typename Expression, bool Vectorizable, TiledEvaluation Tiling>
EIGEN_STRONG_INLINE void TensorExecutor<Expression, GpuDevice, Vectorizable, Tiling>::run(
    const Expression& expr, const GpuDevice& device) {
  TensorEvaluator<Expression, GpuDevice> evaluator(expr, device);
  const bool needs_assign = evaluator.evalSubExprsIfNeeded(nullptr);
  if (needs_assign) {

    const int block_size = device.maxGpuThreadsPerBlock();
    const int max_blocks = device.getNumGpuMultiProcessors() *
                           device.maxGpuThreadsPerMultiProcessor() / block_size;
    const StorageIndex size = array_prod(evaluator.dimensions());
    // Create a least one block to ensure we won't crash when tensorflow calls with tensors of size 0.
    const int num_blocks = numext::maxi<int>(numext::mini<int>(max_blocks, divup<int>(size, block_size)), 1);

    LAUNCH_GPU_KERNEL(
        (EigenMetaKernel<TensorEvaluator<Expression, GpuDevice>, StorageIndex>),
        num_blocks, block_size, 0, device, evaluator, size);
  }
  evaluator.cleanup();
}

#endif  // EIGEN_GPUCC
#endif  // EIGEN_USE_GPU

// SYCL Executor policy
#ifdef EIGEN_USE_SYCL

<<<<<<< HEAD
template <typename Evaluator>
struct ExecExprFunctorKernel {
  typedef typename Evaluator::Index Index;
  Evaluator evaluator;
  const Index range;
  template <typename Scratch>
  EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE ExecExprFunctorKernel(
      const Scratch, Evaluator evaluator_, const Index range_)
      : evaluator(evaluator_), range(range_) {}

  EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE void operator()(
      cl::sycl::nd_item<1> itemID) {
    compute(itemID);
  }
  template <bool is_vec = Evaluator::PacketAccess>
  EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE typename std::enable_if<!is_vec>::type
  compute(const cl::sycl::nd_item<1>& itemID) {
    Index gId = static_cast<Index>(itemID.get_global_linear_id());
    Index total_threads = itemID.get_global_range(0);

=======
template <bool Vectorizable, typename Evaluator>
struct ExecExprFunctorKernel_impl {
  typedef typename Evaluator::Index Index;
  const Index range;
  const Index vectorizable_threads;
  Evaluator evaluator;
  EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE ExecExprFunctorKernel_impl(
      const Index range_, const Index vectorizable_threads_,
      Evaluator evaluator_)
      : range(range_), vectorizable_threads(vectorizable_threads_),
        evaluator(evaluator_) {}

  EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE void
  operator()(cl::sycl::nd_item<1> itemID) {
    Index gId = static_cast<Index>(itemID.get_global_linear_id());
    Index total_threads = itemID.get_global_range(0);
    EIGEN_UNROLL_LOOP
>>>>>>> master
    for (Index i = gId; i < range; i += total_threads) {
      evaluator.evalScalar(i);
    }
  }
<<<<<<< HEAD
  template <bool is_vec = Evaluator::PacketAccess>
  EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE typename std::enable_if<is_vec>::type
  compute(const cl::sycl::nd_item<1>& itemID) {
    const Index vectorizedRange =
        (range / Evaluator::PacketSize) * Evaluator::PacketSize;
    Index gId = static_cast<Index>(itemID.get_global_linear_id());
    const Index step = Evaluator::PacketSize * itemID.get_global_range(0);
    const Index start = Evaluator::PacketSize * gId;
    for (Index i = start; i < vectorizedRange; i += step) {
      evaluator.evalPacket(i);
    }
    gId += vectorizedRange;
    for (Index i = gId; i < range; i += itemID.get_global_range(0)) {
      evaluator.evalScalar(i);
=======
};

template <typename Evaluator>
struct ExecExprFunctorKernel_impl<true, Evaluator> {
  typedef typename Evaluator::Index Index;
  const Index range;
  const Index vectorizable_threads;
  Evaluator evaluator;
  EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE ExecExprFunctorKernel_impl(
      const Index range_, const Index vectorizable_threads_,
      Evaluator evaluator_)
      : range(range_), vectorizable_threads(vectorizable_threads_),
        evaluator(evaluator_) {}

  EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE void
  operator()(cl::sycl::nd_item<1> itemID) {
    Index gId = static_cast<Index>(itemID.get_global_linear_id());
    if (gId < vectorizable_threads) {
      const Index PacketSize = Eigen::internal::unpacket_traits<
          typename Evaluator::PacketReturnType>::size;
      evaluator.evalPacket(gId * PacketSize);
      gId += (vectorizable_threads * PacketSize);
      EIGEN_UNROLL_LOOP
      for (Index i = gId; i < range; i += vectorizable_threads) {
        evaluator.evalScalar(i);
      }
>>>>>>> master
    }
  }
};

<<<<<<< HEAD
template <typename Expression, bool Vectorizable, TiledEvaluation Tiling>
class TensorExecutor<Expression, Eigen::SyclDevice, Vectorizable, Tiling> {
 public:
  typedef typename Expression::Index Index;
  static EIGEN_STRONG_INLINE void run(const Expression& expr,
                                      const Eigen::SyclDevice& dev) {
    typedef Eigen::TensorEvaluator<Expression, Eigen::SyclDevice> Evaluator;
    Evaluator evaluator(expr, dev);
    const bool needs_assign = evaluator.evalSubExprsIfNeeded(NULL);
=======
template <typename Expr, bool NonZeroVectoriseSize, typename Evaluator>
struct ExecExprFunctorKernel
    : ExecExprFunctorKernel_impl<
          ::Eigen::internal::IsVectorizable<Eigen::SyclDevice, Expr>::value,
          Evaluator> {
  ExecExprFunctorKernel(const Index range_, const Index vectorizable_threads_,
                        const Evaluator &evaluator)
      : ExecExprFunctorKernel_impl<
            ::Eigen::internal::IsVectorizable<Eigen::SyclDevice, Expr>::value,
            Evaluator>(range_, vectorizable_threads_, evaluator) {}
};

template <typename Expr, typename Evaluator>
struct ExecExprFunctorKernel<Expr, false, Evaluator>
    : ExecExprFunctorKernel_impl<false, Evaluator> {
  ExecExprFunctorKernel(const Index range_, const Index vectorizable_threads_,
                        const Evaluator &evaluator)
      : ExecExprFunctorKernel_impl<false, Evaluator>(
            range_, vectorizable_threads_, evaluator) {}
};

template <typename Expression, bool Vectorizable, TiledEvaluation Tiling>
class TensorExecutor<Expression, Eigen::SyclDevice, Vectorizable, Tiling> {
  public:
  typedef typename Expression::Index Index;
   static EIGEN_STRONG_INLINE void run(const Expression &expr, const Eigen::SyclDevice &dev) {
    Eigen::TensorEvaluator<Expression, Eigen::SyclDevice> evaluator(expr, dev);
    const bool needs_assign = evaluator.evalSubExprsIfNeeded(nullptr);
>>>>>>> master
    if (needs_assign) {
      Index range, GRange, tileSize;
      Index total_size = ::Eigen::internal::array_prod(evaluator.dimensions());
      total_size = (total_size == 0) ? 1 : total_size;
<<<<<<< HEAD
      const int PacketSize =
          Eigen::PacketType<typename Evaluator::CoeffReturnType,
                            Eigen::SyclDevice>::size;
      Index vectorizable_threads = static_cast<Index>(total_size / PacketSize);
      dev.parallel_for_setup(vectorizable_threads, tileSize, range, GRange);
      range = total_size;

      dev.template nullary_kernel_launcher<
          typename Evaluator::CoeffReturnType,
          ExecExprFunctorKernel<Evaluator> >(
          evaluator,
          cl::sycl::nd_range<1>(cl::sycl::range<1>(GRange),
                                cl::sycl::range<1>(tileSize)),
          Index(1), range);
=======
      const int PacketSize = Eigen::PacketType<
          typename Eigen::TensorEvaluator<Expression, Eigen::SyclDevice>::CoeffReturnType,
          Eigen::SyclDevice>::size;
      Index vectorizable_threads =
          static_cast<Index>(total_size / PacketSize);
      dev.parallel_for_setup(vectorizable_threads, tileSize, range, GRange);
      range = total_size;
      auto f = [&](cl::sycl::handler &cgh) {
        evaluator.bind(cgh);
        typedef ExecExprFunctorKernel<Expression, true,
                                      Eigen::TensorEvaluator<Expression, Eigen::SyclDevice>>
            conditional_vectorized_kernel;

        typedef ExecExprFunctorKernel<Expression, false,
                                      Eigen::TensorEvaluator<Expression, Eigen::SyclDevice>>
            non_vectorized_kernel;
// This is to make sure that an expression with a size less than vectorized size
// will not call the vectorized kernel.
// The reason for having this kernel is that the vectorisable parameter is a
// compile-time parameter,
// however, the size of a tensor is a run-time parameter
        (vectorizable_threads)
            ? cgh.parallel_for(
#ifdef EIGEN_SYCL_USE_PROGRAM_CLASS
                  dev.program().template get_kernel<vectorized_kernel>(),
#endif
                  cl::sycl::nd_range<1>(cl::sycl::range<1>(GRange),
                                        cl::sycl::range<1>(tileSize)),
                  conditional_vectorized_kernel(range, vectorizable_threads,
                                                evaluator))
            : cgh.parallel_for(
#ifdef EIGEN_SYCL_USE_PROGRAM_CLASS
                  dev.program().template get_kernel<non_vectorized_kernel>(),
#endif
                  cl::sycl::nd_range<1>(cl::sycl::range<1>(GRange),
                                        cl::sycl::range<1>(tileSize)),
                  non_vectorized_kernel(range, vectorizable_threads,
                                        evaluator));
      };
      cl::sycl::event e;
      EIGEN_SYCL_TRY_CATCH(e = dev.sycl_queue().submit(f));
      dev.async_synchronize(e);
>>>>>>> master
    }
    evaluator.cleanup();
  }
};

#endif

} // end namespace internal

} // end namespace Eigen

#endif // EIGEN_CXX11_TENSOR_TENSOR_EXECUTOR_H
