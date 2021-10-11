/*
 *  Copyright 2019-2021 Diligent Graphics LLC
 *  Copyright 2015-2019 Egor Yusov
 *  
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *  
 *      http://www.apache.org/licenses/LICENSE-2.0
 *  
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  In no event and under no legal theory, whether in tort (including negligence), 
 *  contract, or otherwise, unless required by applicable law (such as deliberate 
 *  and grossly negligent acts) or agreed to in writing, shall any Contributor be
 *  liable for any damages, including any direct, indirect, special, incidental, 
 *  or consequential damages of any character arising as a result of this License or 
 *  out of the use or inability to use the software (including but not limited to damages 
 *  for loss of goodwill, work stoppage, computer failure or malfunction, or any and 
 *  all other commercial damages or losses), even if such Contributor has been advised 
 *  of the possibility of such damages.
 */

// This file is derived from the open source project provided by Intel Corportaion that
// requires the following notice to be kept:
//--------------------------------------------------------------------------------------
// Copyright 2013 Intel Corporation
// All Rights Reserved
//
// Permission is granted to use, copy, distribute and prepare derivative works of this
// software for any purpose and without fee, provided, that the above copyright notice
// and this statement appear in all copies.  Intel makes no representations about the
// suitability of this software for any purpose.  THIS SOFTWARE IS PROVIDED "AS IS."
// INTEL SPECIFICALLY DISCLAIMS ALL WARRANTIES, EXPRESS OR IMPLIED, AND ALL LIABILITY,
// INCLUDING CONSEQUENTIAL AND OTHER INDIRECT DAMAGES, FOR THE USE OF THIS SOFTWARE,
// INCLUDING LIABILITY FOR INFRINGEMENT OF ANY PROPRIETARY RIGHTS, AND INCLUDING THE
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  Intel does not
// assume any responsibility for any errors which may appear in this software nor any
// responsibility to update it.
//--------------------------------------------------------------------------------------

#pragma once

#include <cassert>
#include <memory>

namespace Diligent
{

// Structure describing quad tree node location
struct QuadTreeNodePosition
{
  // Position in a tree
  int horz_order;
  int vert_order;
  int level;

  QuadTreeNodePosition(int h, int v, int l) : horz_order(h), vert_order(v), level(l)
  {
    VERIFY_EXPR(h < (1 << l)); // h< 2^l
    VERIFY_EXPR(v < (1 << l)); // v< 2^l
  }
  QuadTreeNodePosition() : horz_order(0), vert_order(0), level(0) {}

  // Gets position of a child
  /*
    -----------------
    |       |       |
    |   2   |   3   |
    |       |       |
    -----------------
    |       |       |
    |   0   |   1   |
    |       |       |
    -----------------
    */
  inline friend QuadTreeNodePosition GetChildPosition(const QuadTreeNodePosition& parent, unsigned int birth_order)
  {
    VERIFY_EXPR(birth_order >= 0 && birth_order < 4);
    return QuadTreeNodePosition(parent.horz_order * 2 + (birth_order & 1), parent.vert_order * 2 + (birth_order >> 1), parent.level + 1);
  }

  // Gets location of a parent
  inline friend QuadTreeNodePosition GetParentLocation(const QuadTreeNodePosition& node)
  {
    assert(node.level > 0);
    return QuadTreeNodePosition(node.horz_order / 2, node.vert_order / 2, node.level - 1);
  }
};

// Base class for iterators traversing the quad tree
class HierarchyIteratorBase
{
public:
      operator const QuadTreeNodePosition&() const { return now_; }
  int Level() const { return now_.level; }
  int Horz() const { return now_.horz_order; }
  int Vert() const { return now_.vert_order; }

protected:
  QuadTreeNodePosition now_;
  int                  now_level_size_{1};
};

// Iterator for recursively traversing the quad tree starting from the root up to the specified level
class HierarchyIterator : public HierarchyIteratorBase
{
public:
  HierarchyIterator(int num_levels) : num_levels_(num_levels) {}
  bool IsValid() const { return now_.level < num_levels_; }
  void Next()
  {
    if (++now_.horz_order == now_level_size_)
    {
      now_.horz_order = 0;
      if (++now_.vert_order == now_level_size_)
      {
        now_.vert_order = 0;
        now_level_size_ = 1 << ++now_.level;
      }
    }
  }

private:
  int num_levels_;
};

// Iterator for recursively traversing the quad tree starting from the specified level up to the root
class HierarchyReverseIterator : public HierarchyIteratorBase
{
public:
  HierarchyReverseIterator(int num_levels)
  {
    now_.level      = num_levels - 1;
    now_level_size_ = 1 << now_.level;
  }
  bool IsValid() const { return now_.level >= 0; }
  void Next()
  {
    if (++now_.horz_order == now_level_size_)
    {
      now_.horz_order = 0;
      if (++now_.vert_order == now_level_size_)
      {
        now_.vert_order = 0;
        now_level_size_ = 1 << --now_.level;
      }
    }
  }
};

// Template class for the node of a dynamic quad tree
template <typename NodeDataType>
class DynamicQuadTreeNode
{
public:
  typedef DynamicQuadTreeNode*                 RawPtr;
  typedef std::unique_ptr<DynamicQuadTreeNode> UniPtr;

  DynamicQuadTreeNode() : parent_(NULL) {}

  NodeDataType&       GetData() { return data_; }
  const NodeDataType& GetData() const { return data_; }

  RawPtr GetParent() const { return parent_; }

  void GetChildren(const RawPtr& lb_child, const RawPtr& rb_child, const RawPtr& lt_child, const RawPtr& rt_child) const
  {
    lb_child = lb_child_.get();
    rb_child = rb_child_.get();
    lt_child = lt_child_.get();
    rt_child = rt_child_.get();
  }

  void GetChildren(RawPtr& lb_child, RawPtr& rb_child, RawPtr& lt_child, RawPtr& rt_child)
  {
    lb_child = lb_child_.get();
    rb_child = rb_child_.get();
    lt_child = lt_child_.get();
    rt_child = rt_child_.get();
  }

  // Attahes specified children to the tree
  void CreateChildren(UniPtr lb_child, UniPtr rb_child, UniPtr lt_child, UniPtr rt_child);
  // Creates children UN-ATTACHED to the tree
  void CreateFloatingChildren(UniPtr& lb_child, UniPtr& rb_child, UniPtr& lt_child, UniPtr& rt_child);
  // Destroys ALL children for the node
  void DestroyChildren();

  const QuadTreeNodePosition& GetPos() const { return pos_; }

  void SetPos(const QuadTreeNodePosition& pos) { pos_ = pos; }

private:
  DynamicQuadTreeNode(RawPtr parent, int birth_order) : parent_(parent), pos_(GetChildPosition(parent->pos_, birth_order)) {}

  NodeDataType data_;

  UniPtr lb_child_;
  UniPtr rb_child_;
  UniPtr lt_child_;
  UniPtr rt_child_;
  RawPtr parent_;

  QuadTreeNodePosition pos_;
};

template <typename NodeDataType>
void DynamicQuadTreeNode<NodeDataType>::CreateChildren(UniPtr lb_child, UniPtr rb_child, UniPtr lt_child, UniPtr rt_child)
{
  assert(!lb_child_.get());
  assert(!rb_child_.get());
  assert(!lt_child_.get());
  assert(!rt_child_.get());

  lb_child_ = lb_child;
  rb_child_ = rb_child;
  lt_child_ = lt_child;
  rt_child_ = rt_child;
}

template <typename NodeDataType>
void DynamicQuadTreeNode<NodeDataType>::CreateFloatingChildren(UniPtr& lb_child, UniPtr& rb_child, UniPtr& lt_child, UniPtr& rt_child)
{
  lb_child.reset(new DynamicQuadTreeNode(this, 0));
  rb_child.reset(new DynamicQuadTreeNode(this, 1));
  lt_child.reset(new DynamicQuadTreeNode(this, 2));
  rt_child.reset(new DynamicQuadTreeNode(this, 3));
}

template <typename NodeDataType>
void DynamicQuadTreeNode<NodeDataType>::DestroyChildren()
{
  if (lb_child_.get()) lb_child_->DestroyChildren();
  if (rb_child_.get()) rb_child_->DestroyChildren();
  if (lt_child_.get()) lt_child_->DestroyChildren();
  if (rt_child_.get()) rt_child_->DestroyChildren();

  lb_child_.reset();
  rb_child_.reset();
  lt_child_.reset();
  rt_child_.reset();
}

} // namespace Diligent
