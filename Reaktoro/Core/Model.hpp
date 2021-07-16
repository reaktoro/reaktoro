// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Core/Params.hpp>

namespace Reaktoro {

template<typename Signature>
class Model;

/// The functional signature of functions that evaluates properties.
template<typename ResultRef, typename... Args>
using ModelEvaluator = Fn<void(ResultRef res, Args... args)>;

/// The functional signature of functions that calculates properties.
template<typename Result, typename... Args>
using ModelCalculator = Fn<Result(Args... args)>;

/// The functional signature of functions that serialize a Model object.
using ModelSerializer = Fn<yaml()>;

/// The class used to represent a model function and its parameters.
/// @ingroup Core
template<typename Result, typename... Args>
class Model<Result(Args...)>
{
public:
    /// The reference type of result type.
    /// In case a custom reference type other than `Result&` is needed,
    /// use `REAKTORO_DEFINE_REFERENCE_TYPE_OF(Result, CustomResultRef)`.
    using ResultRef = Ref<Result>;

    /// Construct a default Model function object.
    Model()
    {}

    /// Construct a Model function object with given model evaluator function and its parameters.
    /// @param evalfn The function that evaluates the model.
    /// @param params The parameters of the underlying model function.
    /// @param serializerfn The function that serializes the underlying model function to yaml format.
    Model(const ModelEvaluator<ResultRef, Args...>& evalfn, const Params& params = {}, const ModelSerializer& serializerfn = {})
    : m_params(params), m_evalfn(evalfn), m_serializerfn(serializerfn)
    {
        assert(evalfn);

        m_calcfn = [evalfn](const Args&... args, const Params& w) -> Result
        {
            Result res;
            evalfn(res, args...);
            return res;
        };
    }

    /// Construct a Model function object with given direct model calculator and its parameters.
    /// @param calcfn The function that calculates the model properties and return them.
    /// @param params The parameters of the underlying model function.
    /// @param serializerfn The function that serializes the underlying model function to yaml format.
    Model(const ModelCalculator<Result, Args...>& calcfn, const Params& params = {}, const ModelSerializer& serializerfn = {})
    : m_params(params), m_serializerfn(serializerfn)
    {
        assert(calcfn);

        m_evalfn = [calcfn](ResultRef res, const Args&... args)
        {
            res = calcfn(args...);
        };

        m_calcfn = [calcfn](const Args&... args, const Params& w) -> Result
        {
            return calcfn(args...);
        };
    }

    /// Construct a Model function object with either a model evaluator or a model calculator function.
    /// This constructor exists so that functions that are not wrapped
    /// into an `std::function` object can be used to construct a Model
    /// function object. Without this constructor, an explicit wrap must
    /// be performed by the used. For example,
    /// `Model(ModelCalculator<real(real,real)>([](real T, real P) { return A + B*T + C*T*P; }))`
    /// can be replaced with `Model([](real T, real P) { return A + B*T + C*T*P; })`.
    /// @param f A model evaluator or a model calculator function.
    template<typename Fun, EnableIf<!isFunction<Fun>>...>
    Model(const Fun& f)
    : Model(asFunction(f))
    {}

    /// Return a new Model function object with memoization for the model calculator.
    auto withMemoization() const -> Model
    {
        Model copy = *this;
        copy.m_calcfn = memoizeLast(copy.m_calcfn); // Here, if `m_calcfn` did not consider `const Params&` as argument, memoization would not know when the parameters have been changed externally!
        return copy;
    }

    /// Evaluate the model with given arguments.
    auto apply(ResultRef res, const Args&... args) const -> void
    {
        assert(m_evalfn);
        m_evalfn(res, args...);
    }

    /// Evaluate the model with given arguments and return the result of the evaluation.
    auto operator()(const Args&... args) const -> Result
    {
        assert(m_calcfn);
        return m_calcfn(args..., m_params);
    }

    /// Return true if this Model function object has been initialized.
    auto initialized() const -> bool
    {
        return m_evalfn != nullptr;
    }

    /// Return true if this Model function object has been initialized.
    operator bool() const
    {
        return initialized();
    }

    /// Return the model evaluator function of this Model function object.
    auto evaluatorFn() const -> const ModelEvaluator<ResultRef, Args...>&
    {
        return m_evalfn;
    }

    /// Return the model calculator function of this Model function object.
    auto calculatorFn() const -> const ModelCalculator<Result, Args..., const Params&>&
    {
        return m_calcfn;
    }

    /// Return the function that serializes the underlying model function to yaml format.
    auto serializerFn() const -> const ModelSerializer
    {
        return m_serializerfn;
    }

    /// Return the model parameters of this Model function object.
    auto params() const -> const Params&
    {
        return m_params;
    }

    /// Return serialization of the underlying model function to yaml format.
    auto serialize() const -> yaml
    {
        return m_serializerfn ? m_serializerfn() : yaml{}; // evaluate m_serializerfn because Param objects may have changed
    }

    /// Return a constant Model function object.
    /// @param param The parameter with the constant value always returned by the Model function object.
    static auto Constant(const Param& param) -> Model
    {
        auto calcfn = [param](const Args&... args) { return param; }; // no need to have `const Params&` in the lambda function here. This is added in the constructor call below!
        Params params = { param };
        return Model(calcfn, params);
    }

private:
    /// The parameters used to initialize the underlying model function.
    /// These parameters can be changed externally and affect the model result.
    /// This is possible because their data is wrapped in a shared pointer.
    /// Care must be taken when memoization is applied to the Model object.
    /// Otherwise, the memoized model cannot realize that embedded Param
    /// objects have been externally modified. For example, consider a model
    /// that depends on temperature and pressure only. In this model, one or
    /// more Param objects may have been captured (e.g., via lambda capture).
    /// If temperature and pressure in a new calculation are the same as last
    /// time, but these captured Param objects have been changed externally,
    /// the memoized version of the Model object will return the cached result
    /// (from last calculation). To prevent this, #m_calcfn below have its
    /// functional signature extended with `const Params&`. By doing this, and
    /// passing along #m_params to its call, its memoized version (@see
    /// withMemoization) will be able to detect if these Param objects have
    /// been changed externally.
    Params m_params;

    /// The underlying model function that performs property evaluations.
    ModelEvaluator<ResultRef, Args...> m_evalfn;

    /// The underlying model function that performs property calculations.
    /// Note the added dependency on `const Params&`.
    /// This is needed for proper memoization optimization!
    ModelCalculator<Result, Args..., const Params&> m_calcfn;

    /// The function that serializes the underlying model function to yaml format.
    /// This has to be a function because if we stored the serialization of the
    /// model at construction and the Param objects associated to it changed at
    /// some point later, then the stored serialization would be out of sync
    /// with the Param objects. By storing a function, the serialization can be
    /// computed at any point, say, after the changes in the Param objects.
    ModelSerializer m_serializerfn;
};

/// Return a reaction thermodynamic model resulting from chaining other models.
template<typename Result, typename... Args>
auto chain(const Vec<Model<Result(Args...)>>& models) -> Model<Result(Args...)>
{
    using ResultRef = Ref<Result>;

    const auto evalfns = vectorize(models, RKT_LAMBDA(model, model.evaluatorFn()));
    const auto serializerfns = vectorize(models, RKT_LAMBDA(model, model.serializerFn()));

    auto evalfn = [=](ResultRef res, const Args&... args)
    {
        for(auto i = 0; i < evalfns.size(); ++i)
            evalfns[i](res, args...);
    };

    auto serializerfn = [serializerfns]() -> yaml
    {
        yaml result;
        for(auto i = 0; i < serializerfns.size(); ++i)
            result.push_back(serializerfns[i]());
        return result;
    };

    Params params;
    for(const auto& model : models)
        for(const auto& param : model.params())
            params.append(param);

    return Model<Result(Args...)>(evalfn, params);
}

/// Return a reaction thermodynamic model resulting from chaining other models.
template<typename Signature>
auto chain(const Model<Signature>& model) -> Model<Signature>
{
    return model;
}

/// Return a reaction thermodynamic model resulting from chaining other models.
template<typename Result, typename... Args, typename... Models>
auto chain(const Model<Result(Args...)>& model, const Models&... models) -> Model<Result(Args...)>
{
    Vec<Model<Result(Args...)>> vec = {model, models...};
    return chain(vec);
}

} // namespace Reaktoro
